#include "cgnspanel.h"
#include <QGroupBox>
#include <QGridLayout>
#include <QResizeEvent>
#include <QTextEdit>
#include <QLineEdit>
#include <QLabel>

CgnsPanel::CgnsPanel(QWidget *parent)
    : QWidget{parent}
{
    this->h1 = 3*40;
    this->h2 = 2*40;
    this->h3 = 3*40;
    int y1 = 0;
    int y2 = y1 + this->h1;
    int y3 = y2 + this->h2;
    int y4 = y3 + this->h3;
    int h4 = 300;

    this->nodeDescGroupBox = new QGroupBox(this);
    this->nodeDescGroupBox->setTitle("Node Description");
    this->nodeDescGroupBox->setGeometry( QRect(0,y1,300,h1) );

    this->parentNodeLabel = new QLabel(this);
    this->parentNodeLabel->setText("Parent Node");

    this->parentNodeLineEdit = new QLineEdit(this);

    this->nodeNameLabel = new QLabel(this);
    this->nodeNameLabel->setText("Node Name");

    this->nodeNameLineEdit = new QLineEdit(this);

    this->nodeLabelLabel = new QLabel(this);
    this->nodeLabelLabel->setText("Node Label");

    this->nodeLabelLineEdit = new QLineEdit(this);

    this->nodeDescGridLayout = new QGridLayout(this->nodeDescGroupBox);
    this->nodeDescGridLayout->addWidget(this->parentNodeLabel, 0, 0, 1, 1);
    this->nodeDescGridLayout->addWidget(this->parentNodeLineEdit, 0, 1, 1, 1);
    this->nodeDescGridLayout->addWidget(this->nodeNameLabel, 1, 0, 1, 1);
    this->nodeDescGridLayout->addWidget(this->nodeNameLineEdit, 1, 1, 1, 1);
    this->nodeDescGridLayout->addWidget(this->nodeLabelLabel, 2, 0, 1, 1);
    this->nodeDescGridLayout->addWidget(this->nodeLabelLineEdit, 2, 1, 1, 1);


    this->linkDescGroupBox = new QGroupBox(this);
    this->linkDescGroupBox->setTitle("Link Description");
    this->linkDescGroupBox->setGeometry( QRect(0,y2,300,h2) );

    this->linkFileLabel = new QLabel(this);
    this->linkFileLabel->setText("Link File");

    this->linkFileLineEdit = new QLineEdit(this);

    this->linkNodeLabel = new QLabel(this);
    this->linkNodeLabel->setText("Link Node");

    this->linkNodeLineEdit = new QLineEdit(this);

    this->linkDescGridLayout = new QGridLayout(this->linkDescGroupBox);
    this->linkDescGridLayout->addWidget(this->linkFileLabel, 0, 0, 1, 1);
    this->linkDescGridLayout->addWidget(this->linkFileLineEdit, 0, 1, 1, 1);
    this->linkDescGridLayout->addWidget(this->linkNodeLabel, 1, 0, 1, 1);
    this->linkDescGridLayout->addWidget(this->linkNodeLineEdit, 1, 1, 1, 1);

    this->dataDescGroupBox = new QGroupBox(this);
    this->dataDescGroupBox->setTitle("Data Description");
    this->dataDescGroupBox->setGeometry( QRect(0,y3,300,h3) );

    this->dataTypeLabel = new QLabel(this);
    this->dataTypeLabel->setText("Data Type");
    this->dataTypeLineEdit = new QLineEdit(this);

    this->dimLabel = new QLabel(this);
    this->dimLabel->setText("Dimensions");
    this->dimLineEdit = new QLineEdit(this);

    this->bytesLabel = new QLabel(this);
    this->bytesLabel->setText("Bytes");
    this->bytesLineEdit = new QLineEdit(this);

    this->dataDescGridLayout = new QGridLayout(this->dataDescGroupBox);
    this->dataDescGridLayout->addWidget(this->dataTypeLabel, 0, 0, 1, 1);
    this->dataDescGridLayout->addWidget(this->dataTypeLineEdit, 0, 1, 1, 1);
    this->dataDescGridLayout->addWidget(this->dimLabel, 1, 0, 1, 1);
    this->dataDescGridLayout->addWidget(this->dimLineEdit, 1, 1, 1, 1);
    this->dataDescGridLayout->addWidget(this->bytesLabel, 2, 0, 1, 1);
    this->dataDescGridLayout->addWidget(this->bytesLineEdit, 2, 1, 1, 1);

    this->nodeDataGroupBox = new QGroupBox(this);
    this->nodeDataGroupBox->setTitle("Node Data");
    this->nodeDataGroupBox->setGeometry( QRect(0,y4,300,h4) );

    this->nodeDataTextEdit = new QTextEdit(this);
    //this->nodeDataTextEdit->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOn);

    QVBoxLayout * vbox = new QVBoxLayout(this->nodeDataGroupBox);
    vbox->addWidget(this->nodeDataTextEdit);
    this->nodeDataGroupBox->setLayout(vbox);
}

void CgnsPanel::resizeEvent(QResizeEvent *event)
{
    this->nodeDescGroupBox->resize(event->size().width(),this->h1);
    this->linkDescGroupBox->resize(event->size().width(),this->h2);
    this->dataDescGroupBox->resize(event->size().width(),this->h3);

    this->nodeDataGroupBox->resize(event->size().width(),event->size().height()-this->h1-this->h2-this->h3);
}
