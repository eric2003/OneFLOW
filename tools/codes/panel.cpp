#include "panel.h"
#include "leftpanel.h"
#include "rightpanel.h"
#include "Terminal.h"
#include "explorer.h"
#include "cgnsview.h"
#include "cgnspanel.h"
#include <QSplitter>
#include <QTextEdit>
#include <QVBoxLayout>
#include <QToolButton>
#include <QStackedWidget>
#include <QCoreApplication>
#include <QResizeEvent>
#include <QPushButton>

QList<WidgetPair*> widgetPairList;
bool initWidgetPairListFlag = false;
void InitWidgetPairList()
{
     if ( initWidgetPairListFlag ) return;
     initWidgetPairListFlag = true;
     widgetPairList.append( new ExplorePair() );
     widgetPairList.append( new ProjectPair() );
     widgetPairList.append( new CgnsPair() );
}

QList<WidgetPair*> &GetWidgetPairList()
{
    InitWidgetPairList();

    return widgetPairList;
}

WidgetPair::WidgetPair()
{
}

WidgetPair::~WidgetPair()
{
}

QWidget *WidgetPair::GetLeft()
{
    return nullptr;
}

QWidget *WidgetPair::GetRight()
{
    return nullptr;
}

ExplorePair::ExplorePair()
{
    this->explorer = new Explorer();
    this->splitterV = new QSplitter(Qt::Vertical);

    QTextEdit* pRightTopEdt = new QTextEdit();
    pRightTopEdt->setText(QObject::tr("Right Top Window"));
    this->splitterV->addWidget(pRightTopEdt);
    this->splitterV->addWidget( new Terminal() );

}

ExplorePair::~ExplorePair()
{
    delete this->explorer;
    delete this->splitterV;
}

QWidget *ExplorePair::GetLeft()
{
    return this->explorer;
}

QWidget *ExplorePair::GetRight()
{
    return this->splitterV;
}

ProjectPair::ProjectPair()
{
    this->leftButton = new QPushButton("Project");
    this->rightButton = new QPushButton("btn2");

}

ProjectPair::~ProjectPair()
{
    delete this->leftButton;
    delete this->rightButton;
}

QWidget *ProjectPair::GetLeft()
{
    return this->leftButton;
}

QWidget *ProjectPair::GetRight()
{
    return this->rightButton;
}


CgnsPair::CgnsPair()
{
    this->cgnsView = new CgnsView();
    this->cgnsPanel = new CgnsPanel();
    this->cgnsView->SetCgnsPanel(this->cgnsPanel);
}

CgnsPair::~CgnsPair()
{
    delete this->cgnsView;
    delete this->cgnsPanel;
}

QWidget *CgnsPair::GetLeft()
{
    return this->cgnsView;
}

QWidget *CgnsPair::GetRight()
{
    return this->cgnsPanel;
}

// tmpInitWidgets _tmpInitWidgets;


Panel::Panel(QWidget *parent)
    : QWidget{parent}
{
    this->setupPanelUi();
    //this->setStyleSheet("background-color: rgb(230, 230, 230);");
    this->buttonPanelWidth = 70;

    expToolButton->setText(QCoreApplication::translate("Panel", "Explorer", nullptr));
    prjToolButton->setText(QCoreApplication::translate("Panel", "Project", nullptr));
    gridToolButton->setText(QCoreApplication::translate("Panel", "Grid", nullptr));

    this->expToolButton->setIcon(QIcon(":/search.png"));
    this->expToolButton->setIconSize(QSize(48,48));
    this->expToolButton->setToolButtonStyle(Qt::ToolButtonTextUnderIcon);

    this->prjToolButton->setIcon(QIcon(":/proj.png"));
    this->prjToolButton->setIconSize(QSize(48,48));
    this->prjToolButton->setToolButtonStyle(Qt::ToolButtonTextUnderIcon);

    this->gridToolButton->setIcon(QIcon(":/grid.png"));
    this->gridToolButton->setIconSize(QSize(48,48));
    this->gridToolButton->setToolButtonStyle(Qt::ToolButtonTextUnderIcon);

    this->splitterH = new QSplitter(Qt::Horizontal, this);

    qDebug() << "this->size() = " << this->size();
    qDebug() << "this->splitterH->size() = " << this->splitterH->size();

    this->leftPanel = new LeftPanel();
    this->rightPanel = new RightPanel();

    this->splitterH->addWidget(this->leftPanel);
    this->splitterH->addWidget(this->rightPanel);

    int myWidth = this->width();
    double r1 = 1.0/5.0;
    double r2 = 1.0 - r1;
    int leftWidth = static_cast<int>(myWidth*r1);
    int rightWidth = static_cast<int>(myWidth*r2);

    QList<int> list;
    list.append(leftWidth);
    list.append(rightWidth);
    this->splitterH->setSizes(list);
    qDebug() << "myWidth = " << myWidth << " leftWidth = " << leftWidth << " rightWidth = " << rightWidth;
    qDebug() << "list="<<list;

    qDebug() << "this->splitterH->sizes()=" << this->splitterH->sizes();

    int myCount = this->splitterH->count();
    qDebug() << "myCount=" << myCount;
    for( int i = 0; i < myCount; ++ i )
    {
        QWidget * widget = this->splitterH->widget(i);
        qDebug() << "Panel::Panel widget->size()=" << widget->size() << "widget->windowTitle()="<< widget->windowTitle();
    }

    qDebug() << "this->splitterH->sizes() after addWidget = " << this->splitterH->sizes();

    QObject::connect(this->expToolButton, &QToolButton::clicked, this, &Panel::onExpToolButtonClicked);
    QObject::connect(this->prjToolButton, &QToolButton::clicked, this, &Panel::onPrjToolButtonClicked);
    QObject::connect(this->gridToolButton, &QToolButton::clicked, this, &Panel::onGridToolButtonClicked);
}

void Panel::setupPanelUi()
{
    this->resize(655, 485);
    verticalLayoutWidget = new QWidget(this);
    verticalLayoutWidget->setObjectName("verticalLayoutWidget");
    verticalLayoutWidget->setGeometry(QRect(0, 0, 89, 481));
    verticalLayout = new QVBoxLayout(verticalLayoutWidget);
    verticalLayout->setSpacing(10);
    verticalLayout->setObjectName("verticalLayout");
    verticalLayout->setContentsMargins(5, 10, 5, 5);
    expToolButton = new QToolButton(verticalLayoutWidget);
    expToolButton->setObjectName("expToolButton");

    verticalLayout->addWidget(expToolButton);

    prjToolButton = new QToolButton(verticalLayoutWidget);
    prjToolButton->setObjectName("prjToolButton");

    verticalLayout->addWidget(prjToolButton);

    gridToolButton = new QToolButton(verticalLayoutWidget);
    gridToolButton->setObjectName("gridToolButton");

    verticalLayout->addWidget(gridToolButton);

    verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Policy::Minimum, QSizePolicy::Policy::Expanding);

    verticalLayout->addItem(verticalSpacer);

}

void Panel::resizeEvent(QResizeEvent *event)
{
    qDebug() << "Panel::resizeEvent event->size() = " << event->size();
    int w = event->size().width();
    int h = event->size().height();
    this->verticalLayoutWidget->resize(this->buttonPanelWidth, h);
    this->splitterH->setGeometry( QRect(this->buttonPanelWidth, 0, w-this->buttonPanelWidth-10, h) );
}

void Panel::onExpToolButtonClicked()
{
    this->leftPanel->stack->setCurrentIndex(0);
    this->rightPanel->stack->setCurrentIndex(0);
}

void Panel::onPrjToolButtonClicked()
{
    this->leftPanel->stack->setCurrentIndex(1);
    this->rightPanel->stack->setCurrentIndex(1);
}

void Panel::onGridToolButtonClicked()
{
    this->leftPanel->stack->setCurrentIndex(2);
    this->rightPanel->stack->setCurrentIndex(2);
}

