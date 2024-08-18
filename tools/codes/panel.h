#ifndef PANEL_H
#define PANEL_H

#include <QWidget>

class QSplitter;
class LeftPanel;
class RightPanel;
class QVBoxLayout;
class QSpacerItem;

QList<QWidget*> &GetWidgetListL();
QList<QWidget*> &GetWidgetListR();

class Panel : public QWidget
{
    Q_OBJECT
public:
    explicit Panel(QWidget *parent = nullptr);
protected:
    void resizeEvent(QResizeEvent *event) override;
private:
    void setupPanelUi();
    void onExpToolButtonClicked();
    void onPrjToolButtonClicked();
    void onGridToolButtonClicked();
signals:
private:
    LeftPanel * leftPanel = nullptr;
    RightPanel * rightPanel = nullptr;
private:
    QSplitter *splitterH;
private:
    QWidget *verticalLayoutWidget;
    QVBoxLayout *verticalLayout;
    QToolButton *expToolButton;
    QToolButton *prjToolButton;
    QToolButton *gridToolButton;
    QSpacerItem *verticalSpacer;
    int buttonPanelWidth;
};

#endif // PANEL_H
