#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

class QMenu;
class QAction;
class CfdThread;
class QProcess;
class Terminal;

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

private:
    void closeEvent( QCloseEvent * ev ) override;
protected:
    bool eventFilter(QObject *obj, QEvent *event) override;
    void keyPressEvent(QKeyEvent *event) override;
private:
    void triggerNew();
    void runCFD();
    void runMPI();
    void runTerminal();
private:
    Ui::MainWindow *ui;
    QMenu * menuFile;
    QAction * actNew;
    QAction * actRun;
    QAction * actRunMPI;
    QAction * actTerminal;
private:
    CfdThread *cfdThread;
    QProcess * terminalProcess = nullptr;
    Terminal * terminal = nullptr;
};
#endif // MAINWINDOW_H
