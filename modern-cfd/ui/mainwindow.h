#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class QProcess;
class QTreeWidget;
class QTreeWidgetItem;
class QTextEdit;
class QLineEdit;

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
    void myReadyRead();
    void readOutput();
    void onReadLineData();
    void OnReturn();
    void OnReturnCmd();
    void ExecuteCommand( const QString & cmd );
    void ExecuteCmdCommand( const QString & cmd );
public:
    void displayJsonTree();
    void AnalysisJsonObj( QJsonObject & jsonObj, QTreeWidgetItem * item );
    void AnalysisJsonValue( QJsonValue & value, QTreeWidgetItem * item );
    void SetItem( QTreeWidgetItem * item, QString & keyname, QJsonValue & obj_value );
    void RunCmd();
private:
    Ui::MainWindow *ui;
    QTreeWidget * treeWidget;
    QTextEdit * cmdEdit;
    QLineEdit * lineEdit;
    QLineEdit * lineEditCmd;
    QProcess * process;
    QProcess * procCmd;
};
#endif // MAINWINDOW_H
