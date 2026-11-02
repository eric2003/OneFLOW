#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

class QMenu;
class QAction;

QT_BEGIN_NAMESPACE
namespace Ui {
class MainWindow;
}
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
public:
    void triggerNew();
    void showDialog();
private:
    Ui::MainWindow *ui;
    QMenu * menuFile;
    QAction * actNew;
};
#endif // MAINWINDOW_H
