#include "rightpanel.h"
#include "ui_rightpanel.h"
#include "panel.h"
#include <QStackedWidget>
#include <QPushButton>
#include <QResizeEvent>
#include <QSplitter>
#include <QTextEdit>

RightPanel::RightPanel(QWidget *parent)
    : QWidget(parent)
    , ui(new Ui::RightPanel)
{
    ui->setupUi(this);

    this->stack = new QStackedWidget(this);

    QList<WidgetPair*> wp = ::GetWidgetPairList();

    for( int i = 0; i < wp.count(); ++ i )
    {
        this->stack->addWidget( wp[i]->GetRight() );
    }


    qDebug() << "this->stack->count()="<< this->stack->count();
}

RightPanel::~RightPanel()
{
    delete ui;
}

void RightPanel::resizeEvent(QResizeEvent *event)
{
    qDebug() << "RightPanel::resizeEvent event->size() = " << event->size();
    int width = event->size().width();
    int height = event->size().height();
    this->stack->resize(width, height);
}
