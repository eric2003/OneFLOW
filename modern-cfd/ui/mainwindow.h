#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class QProcess;
class QTreeWidget;
class QTreeWidgetItem;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
public:
    void initMenu();
    void Run();
public:
    void onReadData();
    void readOutput();
    void onReadLineData();
    void displayTree();
    void displayMyTree();
    void displayJsonTree();
    void AnalysisJsonObj( QJsonObject & jsonObj, QTreeWidgetItem * item );
    void AnalysisJsonValue( QJsonValue & value, QTreeWidgetItem * item );
    void SetItem( QTreeWidgetItem * item, QString & keyname, QJsonValue & obj_value );
private:
    Ui::MainWindow *ui;
    QProcess * process;
};
#endif // MAINWINDOW_H
