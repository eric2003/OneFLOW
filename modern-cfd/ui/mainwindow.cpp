/*---------------------------------------------------------------------------*\
OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
Copyright (C) 2017-2022 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
This file is part of OneFLOW.

OneFLOW is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

OneFLOW is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/
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
#include <QTextEdit>
#include <QLineEdit>
#include <iostream>
#include "MyDataBase.h"

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
    this->setWindowTitle( "OneFLOW-CFD Workbench" );
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

    this->lineEdit = new QLineEdit( this );
    this->lineEdit->setStyle( QStyleFactory::create( "windows" ) );
    this->lineEdit->setGeometry( QRect(100, 360, 600, 20) );

    this->lineEditCmd = new QLineEdit( this );
    this->lineEditCmd->setStyle( QStyleFactory::create( "windows" ) );
    this->lineEditCmd->setGeometry( QRect(100, 400, 600, 20) );

    QProcess::connect( this->lineEdit, & QLineEdit::returnPressed, this, & MainWindow::OnReturn );
    QProcess::connect( this->lineEditCmd, & QLineEdit::returnPressed, this, & MainWindow::OnReturnCmd );

    this->cmdEdit = new QTextEdit( this );

    this->cmdEdit->setStyle( QStyleFactory::create( "windows" ) );
    this->cmdEdit->setGeometry( QRect(100, 440, 600, 200) );

    this->treeWidget = new QTreeWidget( this );
    this->treeWidget->setStyle( QStyleFactory::create( "windows" ) );
    this->treeWidget->setGeometry(QRect(100, 50, 600, 300));
    this->treeWidget->setColumnCount( 1 );
    QStringList labels;
    labels << QString::fromLocal8Bit("CFD²ÎÊý");
    this->treeWidget->setHeaderLabels( labels );
    this->treeWidget->header()->setSectionResizeMode(QHeaderView::Stretch);
    this->treeWidget->setContextMenuPolicy(Qt::CustomContextMenu);

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
    //Game Object, property, behavior
    //Componenet Composition
    //Events
    //Message sending and handling
    //Decoupling event sending handling
    //Scene Management
    //Deal with architecture of modern
    //computer with complex
    //combination of CPU and GPU
    //Lighting,Material,Shader
    //Physical-Based Material PBR



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

    this->procCmd = new QProcess( this );
    QObject::connect( this->procCmd, &QProcess::readyRead, this, &MainWindow::myReadyRead);
    QObject::connect( this->procCmd, &QProcess::readyReadStandardOutput, this, &MainWindow::readOutput);

    //"D:/github/OneFLOW/modern-cfd/build/ui/OneFLOW-1D.exe"
    //"D:/github/OneFLOW/modern-cfd/build/bin/OneFLOW-1D.exe"

    QString program = "D:/github/OneFLOW/modern-cfd/build/bin/OneFLOW-1D.exe";
    QStringList args;
    args.append( "cgns_read_modify_cgnsflow" );
    //this->process->start( program, args );
    QStringList args_test;
    args_test << "-uroot" << "-p123456";
    this->process->start( "mysqlsh.exe", args_test );

    this->displayJsonTree();
    this->RunCmd();

    MyDataBase mydb;
    mydb.Run();
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

void MainWindow::RunCmd()
{
    QString myCommand = "dir";

    QStringList argument;
    argument << "/c" << myCommand;
    this->procCmd->start("cmd", QStringList()<<"/c"<<"ping www.baidu.com");
}

void MainWindow::myReadyRead()
{
    qDebug() << Q_FUNC_INFO;
}

void MainWindow::readOutput()
{
    qDebug() << Q_FUNC_INFO;
    //QByteArray qByteRead = this->procCmd->readAllStandardOutput();
    //qDebug() <<  qPrintable( QString::fromLocal8Bit( qByteRead ).trimmed() );
    //if ( qByteRead.contains("enter your name") ) {
    //    this->procCmd->write( QString("myname" + QString("\n")).toLatin1()) ;
    //    qDebug() << this->procCmd->readAllStandardOutput();
    //}

    //this->procCmd->write( QString("myname" + QString("\n") ).toLatin1());
    QByteArray qByteRead = this->procCmd->readAll() + this->procCmd->readAllStandardOutput();
    this->cmdEdit->append( QString::fromLocal8Bit( qByteRead ) );
}

void MainWindow::onReadLineData()
{
    qDebug() << Q_FUNC_INFO;
    while( this->process->canReadLine() )
    {
        QByteArray qByteRead =  this->process->readLine();
        //qDebug() <<  QString::fromLocal8Bit( qByteRead ).trimmed();
        qDebug() <<  qPrintable( QString::fromLocal8Bit( qByteRead ).trimmed() );
    }
}

void MainWindow::OnReturn()
{
    qDebug() << Q_FUNC_INFO;
    QString cmdString = this->lineEdit->text();
    qDebug() << cmdString;
    this->cmdEdit->append( cmdString );
    this->ExecuteCommand( cmdString );
}

void MainWindow::OnReturnCmd()
{
    qDebug() << Q_FUNC_INFO;
    QString cmdString = this->lineEditCmd->text();
    qDebug() << cmdString;
    this->cmdEdit->append( cmdString );
    this->ExecuteCmdCommand( cmdString );
}

void MainWindow::ExecuteCmdCommand( const QString & myCommand )
{
    qDebug() << Q_FUNC_INFO;
    qDebug() << "begin MainWindow::ExecuteCmdCommand myCommand = " << myCommand;
    this->procCmd->start("cmd", QStringList() << "/c" << myCommand );

    this->procCmd->waitForStarted();
    this->procCmd->waitForFinished();
    QByteArray cmd_byte_error = this->procCmd->readAllStandardError();

    QString cmd_error = QString::fromLocal8Bit( cmd_byte_error );
    this->cmdEdit->append( cmd_error );
    qDebug() << "end MainWindow::ExecuteCmdCommand myCommand";
}

void MainWindow::ExecuteCommand( const QString & myCommand )
{
    qDebug() << Q_FUNC_INFO;
    qDebug() << "begin MainWindow::ExecuteCommand myCommand = " << myCommand;

    QStringList strList =  myCommand.split( " " );
    if ( strList.size() == 0 )
    {
        return;
    }
    QStringList args;
    this->procCmd->setProgram( strList[ 0 ] );

    for ( int i = 1; i < strList.size(); ++ i )
    {
        args << strList[ i ];
    }

    this->procCmd->setArguments( args );
    this->procCmd->start();

    this->procCmd->waitForStarted();
    this->procCmd->waitForFinished();
    QByteArray cmd_byte_error = this->procCmd->readAllStandardError();

    QString cmd_error = QString::fromLocal8Bit( cmd_byte_error );
    this->cmdEdit->append( cmd_error );
    qDebug() << "end MainWindow::ExecuteCommand myCommand";
}
