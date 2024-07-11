#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include <QFileDialog>
#include <QProcess>
#include <QKeyEvent>
#include <thread>
#include <iostream>
#include "CfdThread.h"
#include "Terminal.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    this->menuFile = new QMenu("&File", this);
    this->ui->menubar->addMenu( this->menuFile );
    this->actNew = new QAction("New");
    this->actRun = new QAction("Run CFD");
    this->actRunMPI = new QAction("Run MPI");
    this->actTerminal = new QAction("TERMINAL");
    this->menuFile->addAction( this->actNew );
    this->menuFile->addAction( this->actRun );
    this->menuFile->addAction( this->actRunMPI );
    this->menuFile->addAction( this->actTerminal );

    QObject::connect(this->actNew, &QAction::triggered, this, &MainWindow::triggerNew);
    QObject::connect(this->actRun, &QAction::triggered, this, &MainWindow::runCFD);
    QObject::connect( this->actRunMPI, &QAction::triggered, this, &MainWindow::runMPI);
    QObject::connect( this->actTerminal, &QAction::triggered, this, &MainWindow::runTerminal);
    this->cfdThread = new CfdThread();
    //this->terminal = new Terminal(this);
    this->terminal = new Terminal();
    //this->terminal->setWindowTitle("Terminal");
    //this->terminal->setGeometry( QRect(10, 50, 600, 400) );
    this->terminal->show();
    //this->terminal->setFocus();
    qDebug() << "MainWindow::MainWindow";
}

MainWindow::~MainWindow()
{
    delete ui;
    delete this->cfdThread;
    delete this->terminalProcess;
}

void MainWindow::closeEvent( QCloseEvent * ev )
{
    if (this->cfdThread->isRunning())
    {
        this->cfdThread->stopCFD();
        this->cfdThread->wait();
    }
}

bool MainWindow::eventFilter(QObject *obj, QEvent *event)
{
    qDebug() << "MainWindow::eventFilter";
    if (event->type() == QEvent::KeyPress)
    {
        QKeyEvent *keyEvent = static_cast<QKeyEvent *>(event);
        if (keyEvent->key() == Qt::Key_Return)
        {
            qDebug() << "Enter key pressed!";
            // 执行回车键按下后的操作
            return true;
        }
    }

    return QObject::eventFilter(obj, event);
}

void MainWindow::keyPressEvent(QKeyEvent *event)
{
    qDebug() << "MainWindow::keyPressEvent";
}

void MainWindow::triggerNew()
{
    QFileDialog::getOpenFileName(this,"haha",".","*.txt");
}

void MainWindow::runCFD()
{
    this->cfdThread->start();
}

void MainWindow::runMPI()
{
    QProcess process;

    QString program = "C:/Windows/System32/calc.exe";
    QStringList arguments;
    process.start(program, arguments);

    if ( !process.waitForStarted() ) {
        qDebug() << "Failed to start external application";
        return;
    }

    qDebug() << "haha1";

    if ( !process.waitForFinished() ) {
        qDebug() << "External application crashed or timed out";
        return;
    }
    qDebug() << "haha2";
    qDebug() << "Process output:" << process.readAllStandardOutput();
    qDebug() << "External application finished with exit code:" << process.exitCode();
    qDebug() << process.errorString();
}

void MainWindow::runTerminal()
{
    if ( !terminalProcess )
    {
        this->terminalProcess = new QProcess( this );
    }
}



