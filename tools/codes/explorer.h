#ifndef EXPLORER_H
#define EXPLORER_H

#include <QWidget>
class QTreeView;
class QStandardItem;

class Explorer : public QWidget
{
    Q_OBJECT
public:
    explicit Explorer(QWidget *parent = nullptr);
protected:
    void resizeEvent(QResizeEvent *event) override;
private:
    void SetUpExplorerInfo();
    void SetUpHeader();
    void PrintDirInfo();
    void AddExplorerItem(QStandardItem *rootItem);
    void AddExplorerItem(QStandardItem *rootItem, const QString &filename);
signals:
private:
    QTreeView * treeView;
};

#endif // EXPLORER_H
