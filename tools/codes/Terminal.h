#ifndef TERMINAL_H
#define TERMINAL_H

#include <QWidget>
#include <QEvent>
#include <QKeyEvent>
#include <QProcess>

QT_BEGIN_NAMESPACE
namespace Ui {
class Terminal;
}
QT_END_NAMESPACE

class QProcess;
class QTextEdit;

class QTextEdit;

extern QTextEdit *textEdit;

class Terminal : public QWidget
{
    Q_OBJECT
public:
    explicit Terminal(QWidget *parent = nullptr);
protected:
    bool eventFilter(QObject *obj, QEvent *event) override;
    void resizeEvent(QResizeEvent *event) override;
signals:
private:
    Ui::Terminal *ui;
    QProcess * procCmd = nullptr;
private:
    void OnReturnKeyPressed();
    void OnKeyBackPressed();
    void OnKeyUpPressed();
    void ProcessTabKey();

    void Analysis( const QString & cmdString );
    QString GetCommandString();
    void CommandCompletion( const QString &str );
    QString FindCommonString(const QString &stra, const QString &strb);
    QString FindCommonString(const QStringList &strlist);
    bool InStringList(QChar a, int ipos, const QStringList &strlist);
private:
    QString cmdStr;
    QString lastCommand;
    QStringList historyCmdList;
    int word_index=0;
private:
    void OnStarted();
    void ReadOutput();
    void ReadStandardOutput();
    void FinishedProcess();
    void ErrorProcess();
    void OnStateChanged();
};

#endif // TERMINAL_H
