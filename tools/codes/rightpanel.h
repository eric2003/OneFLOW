#ifndef RIGHTPANEL_H
#define RIGHTPANEL_H

#include <QWidget>

class QStackedWidget;
// class QSplitter;
// class Terminal;
// class CgnsPanel;

namespace Ui {
class RightPanel;
}

class RightPanel : public QWidget
{
    Q_OBJECT

public:
    explicit RightPanel(QWidget *parent = nullptr);
    ~RightPanel();

protected:
    void resizeEvent(QResizeEvent *event) override;
public:
    QStackedWidget * stack;
private:
    Ui::RightPanel *ui;
};

#endif // RIGHTPANEL_H
