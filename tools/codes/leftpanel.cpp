#include "leftpanel.h"
#include "ui_leftpanel.h"
#include <QStackedWidget>
#include <QResizeEvent>
#include <QWidget>
#include <QPushButton>
#include "panel.h"

LeftPanel::LeftPanel(QWidget *parent)
    : QWidget(parent)
    , ui(new Ui::LeftPanel)
{
    ui->setupUi(this);

    this->stack = new QStackedWidget(this);
    //this->stack->setStyleSheet("background-color: rgb(243,243,243);");

    QList<WidgetPair*> wp = ::GetWidgetPairList();

    for( int i = 0; i < wp.count(); ++ i )
    {
        this->stack->addWidget( wp[i]->GetLeft() );
    }

    qDebug() << "this->stack->count()="<< this->stack->count();

}

LeftPanel::~LeftPanel()
{
    delete ui;
}

void LeftPanel::resizeEvent(QResizeEvent *event)
{
    // qDebug() << "LeftPanel::resizeEvent event->size() = " << event->size();
    // qDebug() << "LeftPanel::resizeEvent this->stack->size() = " << this->stack->size();
    int width = event->size().width();
    int height = event->size().height();
    this->stack->resize(width, height);
    // qDebug() << "LeftPanel::resizeEvent this->stack->size() = " << this->stack->size();
    // qDebug() << "LeftPanel::resizeEvent width = " << width << " height = " << height;
}
