#ifndef PANEL_H
#define PANEL_H

#include <QWidget>

class QSplitter;
class LeftPanel;
class RightPanel;
class QVBoxLayout;
class QSpacerItem;

class WidgetPair;
QList<WidgetPair*> &GetWidgetPairList();

class WidgetPair
{
public:
    WidgetPair();
    ~WidgetPair();
public:
    virtual QWidget *GetLeft();
    virtual QWidget *GetRight();
};

class Explorer;

class ExplorePair : public WidgetPair
{
public:
    ExplorePair();
    ~ExplorePair();
public:
    QWidget *GetLeft();
    QWidget *GetRight();
public:
    QSplitter *splitterV = nullptr;
    Explorer *explorer = nullptr;
};

class QPushButton;

class ProjectPair : public WidgetPair
{
public:
    ProjectPair();
    ~ProjectPair();
public:
    QWidget *GetLeft();
    QWidget *GetRight();
public:
    QPushButton *leftButton = nullptr;
    QPushButton *rightButton = nullptr;
};

class CgnsView;
class CgnsPanel;
class CgnsPair : public WidgetPair
{
public:
    CgnsPair();
    ~CgnsPair();
public:
    QWidget *GetLeft();
    QWidget *GetRight();
public:
    CgnsView *cgnsView = nullptr;
    CgnsPanel *cgnsPanel = nullptr;
};


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
