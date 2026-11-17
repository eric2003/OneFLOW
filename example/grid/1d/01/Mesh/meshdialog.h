#ifndef MESHDIALOG_H
#define MESHDIALOG_H

#include <QDialog>

namespace Ui {
class MeshDialog;
}

class BC;
class QComboBox;

class MeshDialog : public QDialog
{
    Q_OBJECT

public:
    explicit MeshDialog(QWidget *parent = nullptr);
    ~MeshDialog();

private slots:
    void on_buttonBox_accepted();
private:
    void DumpGrid( const std::string &gridfilename, int nPoints, double xstart, double xend, BC &left, BC &right );
    void AddBcItems(QComboBox *comBox);
private:
    Ui::MeshDialog *ui;
};

#endif // MESHDIALOG_H
