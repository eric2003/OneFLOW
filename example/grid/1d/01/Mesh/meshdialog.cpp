#include "meshdialog.h"
#include "ui_meshdialog.h"
#include "grid.h"

std::map<std::string, BCType_t> bctype_map{
    {"BCTypeNull", BCTypeNull},
    {"BCAxisymmetricWedge", BCAxisymmetricWedge},
    {"BCDegenerateLine", BCDegenerateLine},
    {"BCDegeneratePoint", BCDegeneratePoint},
    {"BCDirichlet", BCDirichlet},
    {"BCExtrapolate", BCExtrapolate},
    {"BCFarfield", BCFarfield},
    {"BCGeneral", BCGeneral},
    {"BCInflow", BCInflow},
    {"BCInflowSubsonic", BCInflowSubsonic},
    {"BCOutflowSupersonic", BCOutflowSupersonic},
    {"BCSymmetryPlane", BCSymmetryPlane},
    {"BCSymmetryPolar", BCSymmetryPolar},
    {"BCTunnelInflow", BCTunnelInflow},
    {"BCTunnelOutflow", BCTunnelOutflow},
    {"BCWall", BCWall},
    {"BCWallInviscid", BCWallInviscid},
    {"BCWallViscous", BCWallViscous},
    {"BCWallViscousHeatFlux", BCWallViscousHeatFlux},
    {"BCWallViscousIsothermal", BCWallViscousIsothermal},
    {"FamilySpecified", FamilySpecified}
};

MeshDialog::MeshDialog(QWidget *parent)
    : QDialog(parent)
    , ui(new Ui::MeshDialog)
{
    ui->setupUi(this);
    this->AddBcItems(ui->comboBoxBcL);
    this->AddBcItems(ui->comboBoxBcR);
}

MeshDialog::~MeshDialog()
{
    delete ui;
}

void MeshDialog::AddBcItems(QComboBox *comBox)
{
    for( std::map<std::string, BCType_t>::iterator iter = bctype_map.begin(); iter != bctype_map.end(); ++ iter )
    {
        comBox->addItem(QString::fromStdString(iter->first));
    }

}

void MeshDialog::on_buttonBox_accepted()
{
    qDebug() << "MeshDialog::on_buttonBox_accepted()";
    QString txt = ui->npoints_lineEdit->text();
    int nPoints = txt.toInt();
    qDebug() << "nPoints = " << nPoints;
    double xstart = ui->lineEdit_Start_X->text().toDouble();
    double xend = ui->lineEdit_End_X->text().toDouble();

    std::string gridfilename = ui->lineEditGrid->text().toStdString();
    BC left, right;
    left.bcName = ui->lineEditBcNameL->text().toStdString();
    std::map<std::string, BCType_t>::iterator iterL = bctype_map.find(ui->comboBoxBcL->currentText().toStdString());
    left.bctype = iterL->second;
    right.bcName = ui->lineEditBcNameR->text().toStdString();
    std::map<std::string, BCType_t>::iterator iterR = bctype_map.find(ui->comboBoxBcR->currentText().toStdString());
    right.bctype = iterR->second;

    DumpGrid( gridfilename, nPoints, xstart, xend, left, right );
}

void MeshDialog::DumpGrid( const std::string &gridfilename, int nPoints, double xstart, double xend, BC &left, BC &right )
{
    write_grid_str( gridfilename, nPoints, xstart, xend );
    write_bc_str( gridfilename, left, right );
}

