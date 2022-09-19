#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include <QMenuBar>
#include <QMessageBox>
#include <QProcess>
#include <QDir>
#include <QTreeWidget>
#include <QHeaderView>
#include <QStyleFactory>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <iostream>

void PrintDateInfo()
{
    const char *date = __DATE__;
    const char *table[] = {
        "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
    };
    int month = std::find(table, table + 12, std::string(date, 3)) - table + 1;
    int day = std::stoi(std::string(date + 4, 2));
    int year = std::stoi(std::string(date + 7, 4));
    //Sep 15 2022
    std::string str = std::string( date, 3 );
    std::string mydate =  __DATE__;
    std::cout << " mydate = " << mydate << "\n";
    std::cout << " date = " << date << "\n";
    std::cout << " str = " << str << "\n";
    std::cout << " month = " << month << "\n";
    std::cout << " day = " << day << "\n";
    std::cout << " year = " << year << "\n";
}

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    this->process = 0;
    this->initMenu();
    this->setWindowTitle( "OneFLOW-CFD GUI" );
    this->Run();
}

MainWindow::~MainWindow()
{
    this->process->close();
    delete this->process;
    delete ui;
}

void MainWindow::initMenu()
{
    QMenu *pFile = new QMenu( tr("File") );
    QMenu *pEdit = new QMenu( tr("Edit") );
    QMenu *pView = new QMenu( tr("View") );
    QMenuBar * pMenuBar = this->menuBar();
    pMenuBar->addMenu( pFile );
    pMenuBar->addMenu( pEdit );
    pMenuBar->addMenu( pView );

    QAction *pAction = new QAction( tr("New"), pFile );
    pFile->addAction( pAction );
    pAction->setCheckable(false);
    pAction->setShortcut( QKeySequence( ("Ctrl+N") ) );

    pAction = new QAction( tr("Open"), pFile );
    pFile->addAction( pAction );

    pAction->setCheckable(false);
    pAction->setShortcut( QKeySequence( ("Ctrl+O") ) );

    pAction = new QAction( tr("Close"), pFile );
    pFile->addAction( pAction );
    //Tool Layer
    //Function Layer
    //Resource Layer
    //Core Layer
    //Platform Layer
    //Scene and Level
    //Script and Graph
    //Game Logic Data
    //Animation,Physics,Render,Script,Camera
    //Operation Systems
    //Platform File Systems
    //Graphics API
    //Platform SDK
    //Resource Importing
    //Runtime Asset Manager
    //Runtime Resource Management
    //A virtual file system to load/unload assets by path reference
    //Manage asset lifespan and reference by handle system
    //Mesh,Material,Texture,Skeleton,Animation
    //internal,external
    //tickLogic, tickRender
    //Camera,Motor,Controller,Animation,Physics
    //Core-Memory Management
    //Memory Pool / Allocator
    //Reduce cache miss
    //Memory alignment
    //Polymorphic Memory Resource(PMR)
    //Cache locality/diffusion
    //Memory Arena
    //Platform-Graphics API
    //Render Hardware Interface(RHI)
    //Transparent different GPU architectures abd SDK
    //Automatic optimization of target platforms
    //Unleash the Creativity
    //Flexible of coding languages


}

void MainWindow::Run()
{
    ::PrintDateInfo();

    std::string str = QCoreApplication::applicationDirPath().toStdString();
    qDebug() << " str = " << str.data() << "\n";

    QDir dir;
    dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks);
    dir.setSorting(QDir::Size | QDir::Reversed);

    QFileInfoList list = dir.entryInfoList();
    std::cout << "     Bytes Filename" << std::endl;
    for ( int i = 0; i < list.size(); ++ i )
    {
        QFileInfo fileInfo = list.at(i);
        std::cout << qPrintable(QString("%1 %2").arg(fileInfo.size(), 10)
            .arg(fileInfo.fileName()));
        std::cout << std::endl;
    }
    std::string cfd_exe = "OneFLOW-1D.exe";
    QString abs_cfd_exe = dir.absoluteFilePath( QString::fromStdString(cfd_exe) );
    std::cout << "abs_cfd_exe = " << abs_cfd_exe.toStdString() << std::endl;

    this->process = new QProcess( this );
    QObject::connect( this->process, &QProcess::readyReadStandardOutput, this, &MainWindow::onReadLineData );
    //"D:/github/OneFLOW/modern-cfd/build/ui/OneFLOW-1D.exe"
    //"D:/github/OneFLOW/modern-cfd/build/bin/OneFLOW-1D.exe"

    QString program = "D:/github/OneFLOW/modern-cfd/build/bin/OneFLOW-1D.exe";
    QStringList args;
    args.append( "cgns_read_modify_cgnsflow" );
    this->process->start( program, args );
    //this->displayTree();
    //this->displayMyTree();
    this->displayJsonTree();
    
}

void MainWindow::onReadData()
{
    qDebug() << "MainWindow::onReadData()\n";
    qDebug() << this->process->readAllStandardOutput();
}

void MainWindow::readOutput()
{
    qDebug() << "MainWindow::readOutput()\n";
    while( this->process->canReadLine() )
    {
        qDebug() << this->process->readLine() << "\n";
    }
}

void MainWindow::onReadLineData()
{
    qDebug() << "MainWindow::onReadData()\n";
    while( this->process->canReadLine() )
    {
        QByteArray qByteRead =  this->process->readLine();
        //qDebug() <<  QString::fromLocal8Bit( qByteRead ).trimmed();
        qDebug() <<  qPrintable( QString::fromLocal8Bit( qByteRead ).trimmed() );
    }
}

void MainWindow::displayTree()
{
    QTreeWidget *treeWidget = new QTreeWidget(this);
    treeWidget->setGeometry(QRect(100, 100, 300, 300));
    //treeWidget->setColumnCount(1);
    treeWidget->setColumnCount(2);
    QList<QTreeWidgetItem *> items;
    for ( int i = 0; i < 10; ++ i )
    {
        items.append(new QTreeWidgetItem(static_cast<QTreeWidget *>(nullptr), QStringList(QString("item: %1").arg(i))));
    }

    treeWidget->insertTopLevelItems(0, items);
    //treeWidget->insertTopLevelItems(1, items);

    QList<QTreeWidgetItem *> items1;

    for ( int i = 0; i < 10; ++ i )
    {
        items1.append(new QTreeWidgetItem(static_cast<QTreeWidget *>(nullptr), QStringList(QString("item: %1").arg(100+i))));
    }
        
    //treeWidget->insertTopLevelItems(1, items1);
}

void MainWindow::displayMyTree()
{
    QTreeWidget *treeWidget = new QTreeWidget(this);
    treeWidget->setStyle( QStyleFactory::create( "windows" ) );
    treeWidget->setGeometry(QRect(100, 100, 600, 300));
    treeWidget->setColumnCount( 2 );
    treeWidget->setColumnWidth(0, 200 );
    treeWidget->setColumnWidth(1, 200 );
    QStringList labels;
    labels << "FileName" << "FilePath";
    treeWidget->setHeaderLabels( labels );
    treeWidget->header()->setSectionResizeMode(QHeaderView::Stretch);
    treeWidget->setContextMenuPolicy(Qt::CustomContextMenu);

    QStyle *style = treeWidget->style();
    QIcon bookmarkIcon;
    bookmarkIcon.addPixmap( style->standardPixmap( QStyle::SP_FileIcon ) );

    QTreeWidgetItem * root = new QTreeWidgetItem( treeWidget );
    treeWidget->addTopLevelItem( root );
    root->setIcon( 0, bookmarkIcon );
    root->setExpanded( true );
    root->setText( 0, "RootFile" );
    root->setText( 1, "Root File Path" );

    QTreeWidgetItem * child1 = new QTreeWidgetItem();
    root->addChild( child1 );
    child1->setText( 0, "Child1File" );
    child1->setText( 1, "Child1 File Path" );

    QTreeWidgetItem * child2 = new QTreeWidgetItem();
    root->addChild( child2 );
    child2->setText( 0, "Child2File" );
    child2->setText( 1, "Child2 File Path" );

    QTreeWidgetItem * root1 = new QTreeWidgetItem( treeWidget );
    treeWidget->addTopLevelItem( root1 );
    root1->setExpanded( true );
    root1->setText( 0, "Root1File" );
    root1->setText( 1, "Root1 File Path" );

    QTreeWidgetItem * child11 = new QTreeWidgetItem();
    root1->addChild( child11 );
    child11->setText( 0, "Child11File" );
    child11->setText( 1, "Child11 File Path" );

    QTreeWidgetItem * child12 = new QTreeWidgetItem();
    root1->addChild( child12 );
    child12->setText( 0, "Child12File" );
    child12->setText( 1, "Child12 File Path" );

}

void AnalysisJsonValue( QJsonValue & value );
void AnalysisJsonObj( QJsonObject & jsonObj );

void AnalysisJsonObj( QJsonObject & jsonObj )
{
    std::cout << "\n{\n";
    for ( int i = 0; i < jsonObj.size(); ++ i )
    {
        QString keyname = jsonObj.keys()[i];
        std::cout << keyname.toStdString() << " : ";
        QJsonValue obj_value = jsonObj.value( keyname );
        AnalysisJsonValue( obj_value );
    }
    std::cout << "\n}";
}

void AnalysisJsonValue( QJsonValue & value )
{
    if ( value.isArray() )
    {
        QJsonArray arr = value.toArray();
        std::cout << "\n[\n";
        for ( int i = 0; i < arr.size(); ++ i )
        {
            QJsonValue jv = arr[ i ];
            AnalysisJsonValue( jv );
        }
        std::cout << "\n]\n";
    }
    else if ( value.isObject() )
    {
        QJsonObject jsonObj = value.toObject();
        AnalysisJsonObj( jsonObj );
    }
    else if ( value.isString() )
    {
        QString value_name = value.toString();
        std::cout << value_name.toStdString() << "\n";
    }
    else if ( value.isDouble() )
    {
        double value_value = value.toDouble();
        std::cout << value_value << "\n";
    }
}

void MainWindow::displayJsonTree()
{
    QString jsonFile = "D:/github/OneFLOW/modern-cfd/project_tmp/cgns_read_modify_cgnsflow/script/cfd.json";
    QFile myFile( jsonFile );
    if ( !myFile.open( QFile::ReadOnly ) )
    {
        qWarning("couldn't open json file.");
    }
    QByteArray data = myFile.readAll();
    QJsonParseError err;
    QJsonDocument json_doc = QJsonDocument::fromJson( data, & err );
    QTreeWidget *treeWidget = new QTreeWidget(this);
    if ( json_doc.isNull() )
    {
        qDebug() << err.errorString();
    }
    QJsonObject rootObj = json_doc.object();
    int nLen = rootObj.length();
    qDebug() << " nLen = " << nLen << "\n";
    qDebug() << " rootObj.keys() = " << rootObj.keys() << "\n";
    QList<QString> mykeys = rootObj.keys();

    QByteArray request_body = json_doc.toJson( QJsonDocument::JsonFormat::Indented );
    qDebug() << request_body;

    treeWidget->setStyle( QStyleFactory::create( "windows" ) );
    treeWidget->setGeometry(QRect(100, 100, 600, 300));
    treeWidget->setColumnCount( 1 );
    QStringList labels;
    labels << QString::fromLocal8Bit("CFD²ÎÊý");
    treeWidget->setHeaderLabels( labels );
    treeWidget->header()->setSectionResizeMode(QHeaderView::Stretch);
    treeWidget->setContextMenuPolicy(Qt::CustomContextMenu);

    QTreeWidgetItem * current_item = new QTreeWidgetItem();

    treeWidget->addTopLevelItem( current_item );
    current_item->setExpanded( true );
    current_item->setText( 0, "JSON" );

    this->AnalysisJsonObj( rootObj, current_item );
}

void MainWindow::AnalysisJsonObj( QJsonObject & jsonObj, QTreeWidgetItem * root )
{
    QList<QString> mykeys = jsonObj.keys();
    for ( int i = 0; i < mykeys.length(); ++ i )
    {
        QString keyname = mykeys[ i ];

        QTreeWidgetItem * current_item = new QTreeWidgetItem();
        if ( root )
        {
            root->addChild( current_item );
        }
        
        QJsonValue obj_value = jsonObj.value( keyname );
 
        this->SetItem( current_item, keyname, obj_value );
        this->AnalysisJsonValue( obj_value, current_item );

    }
}

void MainWindow::SetItem( QTreeWidgetItem * item, QString & keyname, QJsonValue & value )
{
    QString text = keyname;
    if ( value.isObject() )
    {
        text.prepend( "{ } " );
    }
    else if ( value.isArray() )
    {
        text.prepend( "[ ] " );
    }
    else if ( value.isString() )
    {
        QString value_str = value.toString();
        text.append( " : " );
        text.append( "\"" );
        text.append( value_str );
        text.append( "\"" );
    }
    else if ( value.isDouble() )
    {
        double double_value = value.toDouble();
        text.append( " : " );
        text.append( QString::number(double_value) );
    }

    item->setExpanded( true );
    item->setText( 0, text );
}

void MainWindow::AnalysisJsonValue( QJsonValue & value, QTreeWidgetItem * root )
{
    if ( value.isObject() )
    {
        QJsonObject jsonObj = value.toObject();
        this->AnalysisJsonObj( jsonObj, root );
    }
    else if ( value.isArray() )
    {
        QJsonArray arr = value.toArray();
        for ( int i = 0; i < arr.size(); ++ i )
        {
            QTreeWidgetItem * current_item = new QTreeWidgetItem();
            if ( root )
            {
                root->addChild( current_item );
            }
            QJsonValue jv = arr[ i ];
            QString keyname = QString::number(i);

            this->SetItem( current_item, keyname, jv );
            AnalysisJsonValue( jv, root );
        }
    }
    else if ( value.isString() )
    {
    }
    else if ( value.isDouble() )
    {
    }
}