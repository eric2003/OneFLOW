#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include <QFileDialog>
#include <QProcess>
#include <QKeyEvent>
#include <QSplitter>
#include <QTextEdit>
#include <QTreeView>
#include <QStandardItemModel>
#include <QHeaderView>
#include "CfdThread.h"
#include "Terminal.h"
#include "panel.h"

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
    //this->ui->statusbar->setStyleSheet("background-color: rgb(0, 255, 0);");
    this->ui->statusbar->setStyleSheet("background-color: rgb(0, 122, 204);");
    qDebug() << "MainWindow::MainWindow";

    this->panel = new Panel(this);
    //this->panel = new Panel();
    this->panel->setGeometry( QRect(0, 25, 600, 400) );
    //this->panel->setStyleSheet("background-color: rgb(240, 240, 240);");
    //this->panel->show();
}

MainWindow::~MainWindow()
{
    delete ui;
    delete this->cfdThread;
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

void MainWindow::resizeEvent(QResizeEvent *event)
{
    qDebug() << "MainWindow::resizeEvent event->size() = " << event->size();
    qDebug() << "MainWindow::resizeEvent this->ui->statusbar->size()=" << this->ui->statusbar->size();
    qDebug() << "MainWindow::resizeEvent this->ui->statusbar->pos()=" << this->ui->statusbar->pos();

    int ypos = this->ui->statusbar->pos().y();

    int splitterTop = 25;
    int splitterHeight = ypos - splitterTop;
    this->panel->setGeometry( QRect(0, splitterTop, this->width(), splitterHeight) );
}

void MainWindow::triggerNew()
{
    //QFileDialog::getOpenFileName(this,"haha",".","*.txt");
    int width = this->size().width();
    int height = this->size().height();
    qDebug() << "MainWindow::triggerNew() this->size()=" << this->size();
    this->ui->statusbar->showMessage("hello",5000);
    qDebug() << "MainWindow::triggerNew() this->ui->statusbar->size()=" << this->ui->statusbar->size();
    qDebug() << "MainWindow::triggerNew() this->ui->statusbar->pos()=" << this->ui->statusbar->pos();

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
}



