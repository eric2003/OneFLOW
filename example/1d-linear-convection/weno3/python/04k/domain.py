# domain.py
from mesh import Mesh

class Domain:
    """计算域：管理物理区域、ghost层、索引映射等逻辑，依赖 Mesh 提供几何信息"""
    def __init__(self, config, mesh):
        """
        初始化计算域
        :param mesh: Mesh实例（静态网格属性）
        :param config: CfdConfig实例（包含recon_scheme/spatial_order）
        """
        self.config = config
        self.mesh = mesh        
        
        # 核心：根据重建格式动态计算nghosts
        self.nghosts = self._calc_nghosts()
        
        # 基于nghosts推导索引
        self.ist = self.nghosts  # 物理网格起始索引
        self.ied = self.ist + mesh.ncells  # 物理网格结束索引
        self.ntcells = mesh.ncells + 2 * self.nghosts  # 总网格数（含ghost）
       
        # 可选：调试信息（可后续移除）
        # print(f"mesh.ncells={mesh.ncells}")
        # print(f"self.config.spatial_order={self.config.spatial_order}")
        # print(f"self.nghosts={self.nghosts}")
        # print(f"self.ist={self.ist}")
        # print(f"self.ied={self.ied}")
        
    def _calc_nghosts(self):
        """内部方法：根据重建格式和阶数计算ghost层数量"""
        scheme = self.config.recon_scheme
        order = self.config.spatial_order
        
        if scheme is None:
            raise ValueError("必须先通过 with_reconstruction 设置重建格式！")
        
        if scheme == "eno":
            nghosts = order
        elif scheme == "weno":
            nghosts = order // 2 + 1
        else:
            raise ValueError(f"未知重建格式 {scheme}，无法计算ghost层！")
        
        if nghosts <= 0:
            raise ValueError(f"计算得到的ghost层数量无效：{nghosts}（阶数{order}，格式{scheme}）")
        return nghosts        
        
    def is_physical_cell(self, idx):
        """判断索引是否在物理网格范围内"""
        return self.ist <= idx < self.ied
    
    def get_physical_indices(self):
        """返回物理网格的索引范围（可直接用于循环）"""
        return range(self.ist, self.ied)