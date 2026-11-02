#include "mainwindow.h"
#include <QApplication>
#include <QStyleFactory>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    a.setStyle(QStyleFactory::create("windows"));
    MainWindow w;
    w.setWindowTitle( "OneFLOW-CFD" );
    w.showMaximized();
    return a.exec();
}
