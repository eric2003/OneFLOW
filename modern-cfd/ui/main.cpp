#include "mainwindow.h"
#include <QApplication>
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

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    ::PrintDateInfo();

    std::string str = QCoreApplication::applicationDirPath().toStdString();
    qDebug() << " str = " << str.data() << "\n";

    MainWindow w;
    w.show();
    return a.exec();
}
