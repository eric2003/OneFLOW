/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "GridPara.h"
#include "DataBase.h"
#include <string>

BeginNameSpace( ONEFLOW )

GridPara grid_para;

void GridPara::Init()
{
    // Prefer the single typed loader, then mirror into legacy fields so
    // existing call sites keep working without a big-bang rewrite.
    const GridConfig cfg = GridConfig::FromDataBase();

    this->gridFile        = cfg.sourceFile;
    this->bcFile          = cfg.bcFile;
    this->targetFile      = cfg.targetFile;
    this->filetype        = std::string( ToString( cfg.sourceType ) );
    this->target_filetype = std::string( ToString( cfg.targetType ) );
    this->topo            = cfg.topo;
    this->multiBlock      = cfg.multiBlock;
    this->gridObj         = static_cast< int >( cfg.objective );
    this->gridScale       = cfg.scale;
    this->axis_dir        = cfg.axisDir;

    this->gridTrans.resize( 3 );
    for ( size_t i = 0; i < 3; ++i )
    {
        this->gridTrans[ i ] = cfg.translate[ i ];
    }
}

GridConfig GridPara::ToConfig() const
{
    GridConfig cfg;
    cfg.objective  = this->objective();
    cfg.sourceType = this->sourceType();
    cfg.targetType = this->targetType();
    cfg.sourceFile = this->gridFile;
    cfg.bcFile      = this->bcFile;
    cfg.targetFile  = this->targetFile;
    cfg.topo        = this->topo;
    cfg.multiBlock  = this->multiBlock;
    cfg.axisDir     = this->axis_dir;
    cfg.scale       = this->gridScale;
    cfg.translate   = { 0.0, 0.0, 0.0 };
    for ( size_t i = 0; i < 3 && i < this->gridTrans.size(); ++i )
    {
        cfg.translate[ i ] = this->gridTrans[ i ];
    }
    return cfg;
}

int GetGridTopoType()
{
    return 0;
}

EndNameSpace
