#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include <QMenuBar>
#include <QMessageBox>
#include <QProcess>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    this->process = 0;
    this->initMenu();
    this->setWindowTitle( "OneFLOW-CFD GUI" );
}

MainWindow::~MainWindow()
{
    this->process->close();
    delete this->process;
    delete ui;
}

void MainWindow::initMenu()
{
    //QProcess::execute( "notepad.exe" );
    //QProcess::execute( "calc.exe" );
    this->process = new QProcess( this );
    process->start( "calc.exe" );
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

 /*   this->menuBar()->addMenu(QString::fromLocal8Bit("文件")); 
    this->menuBar()->addMenu(QString::fromLocal8Bit("编辑"));

    qDebug() << "this->ui->menubar=" << this->ui->menubar << "\n";
    QAction * myAc1 = new QAction(this);  
    myAc1->setText(QString::fromLocal8Bit("新建"));
    myAc1->setStatusTip("This is ac1."); 
    qDebug() << "this->ui->menubar=" << this->ui->menubar << "\n";
    this->ui->menubar->addAction(myAc1); 

    QMessageBox msgBox;
    msgBox.setText("The document has been modified.");
    msgBox.exec();*/

    //QMenuBar * pMenuBar = new QMenuBar( this );
    //if ( !pMenuBar )
    //{
    //    return;
    //}
    //this->setMenuBar( pMenuBar );
    //    
    //QMenu *pFile = new QMenu( tr("File") );
    //pMenuBar->addMenu(pFile);
}

