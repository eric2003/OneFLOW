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

class MyEventFilter : public QObject
{
    Q_OBJECT

protected:
    bool eventFilter(QObject *obj, QEvent *event) override
    {
        //qDebug() << "MyEventFilter::eventFilter";
        if (event->type() == QEvent::KeyPress)
        {
            QKeyEvent *keyEvent = static_cast<QKeyEvent *>(event);
            if (keyEvent->key() == Qt::Key_Return)
            {
                qDebug() << "MyEventFilter Enter key pressed!";
                // 执行回车键按下后的操作
                return true;
            }
        }

        return QObject::eventFilter(obj, event);
    }
};


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
    void ReadOutput();
    void ReadStandardOutput();
    void FinishedProcess();
    void ErrorProcess();
    void OnReturnKeyPressed();
    void OnStateChanged();
    void OnKeyBackPressed();
    void ProcessTabKey();
    QString GetCommandString();
    void CommandCompletion( const QString &str );
private:
    QString FindCommonString(const QString &stra, const QString &strb);
    QString FindCommonString(const QStringList &strlist);
    bool InStringList(QChar a, int ipos, const QStringList &strlist);
private:
    void Analysis( QString & cmdString );
private:
    Ui::Terminal *ui;
    QProcess * procCmd = nullptr;
    QString currentPath;
    QString lastCommand;
    QStringList wordList;
    int word_index=0;
};

#endif // TERMINAL_H
