#ifndef CGNSPANEL_H
#define CGNSPANEL_H

#include <QWidget>
class QGroupBox;
class QTextEdit;
class QGridLayout;
class QLabel;
class QLineEdit;

class CgnsPanel : public QWidget
{
    Q_OBJECT
public:
    explicit CgnsPanel(QWidget *parent = nullptr);
protected:
    void resizeEvent(QResizeEvent *event) override;

signals:
private:
    QGroupBox *nodeDescGroupBox;
    QGroupBox *linkDescGroupBox;
    QGroupBox *dataDescGroupBox;
    QGroupBox *nodeDataGroupBox;
    QGridLayout * nodeDescGridLayout;
    QLabel *parentNodeLabel;
    QLineEdit *parentNodeLineEdit;

    QLabel *nodeNameLabel;
    QLineEdit *nodeNameLineEdit;

    QLabel *nodeLabelLabel;
    QLineEdit *nodeLabelLineEdit;

    QGridLayout * linkDescGridLayout;

    QLabel *linkFileLabel;
    QLineEdit *linkFileLineEdit;

    QLabel *linkNodeLabel;
    QLineEdit *linkNodeLineEdit;

    QGridLayout * dataDescGridLayout;
    QLabel *dataTypeLabel;
    QLineEdit *dataTypeLineEdit;

    QLabel *dimLabel;
    QLineEdit *dimLineEdit;

    QLabel *bytesLabel;
    QLineEdit *bytesLineEdit;

    QTextEdit *nodeDataTextEdit;
private:
    int h1;
    int h2;
    int h3;
};

#endif // CGNSPANEL_H
