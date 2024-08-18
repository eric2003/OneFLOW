#ifndef SIDEPANEL_H
#define SIDEPANEL_H

#include <QWidget>

class QStackedWidget;

namespace Ui {
class LeftPanel;
}

class LeftPanel : public QWidget
{
    Q_OBJECT

public:
    explicit LeftPanel(QWidget *parent = nullptr);
    ~LeftPanel();

protected:
    void resizeEvent(QResizeEvent *event) override;

private slots:

public:
    QStackedWidget * stack;
private:
    Ui::LeftPanel *ui;
};

#endif // SIDEPANEL_H
