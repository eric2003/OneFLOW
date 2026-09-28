#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include "meshdialog.h"


MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    this->menuFile = new QMenu("&File", this);
    this->ui->menubar->addMenu( this->menuFile );
    this->actNew = new QAction("New");
    this->menuFile->addAction( this->actNew );
    QObject::connect(this->actNew, &QAction::triggered, this, &MainWindow::triggerNew);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::triggerNew()
{
    showDialog();
}

void MainWindow::showDialog()
{
    MeshDialog dialog(this);
    dialog.exec();  // 显示对话框，阻塞直到关闭
}
