import sympy as sp

# 定义符号
i = sp.symbols('i', integer=True)
v = sp.Function('v')

# 输入：已知的展开式（不给出原始表达式）
beta0_expanded = (10*v(i)**2/3 - 31*v(i)*v(i+1)/3 + 11*v(i)*v(i+2)/3 + 
                  25*v(i+1)**2/3 - 19*v(i+1)*v(i+2)/3 + 4*v(i+2)**2/3)

beta1_expanded = (13*v(i)**2/3 - 13*v(i)*v(i-1)/3 - 13*v(i)*v(i+1)/3 + 
                  4*v(i-1)**2/3 + 5*v(i-1)*v(i+1)/3 + 4*v(i+1)**2/3)

beta2_expanded = (10*v(i)**2/3 + 11*v(i)*v(i-2)/3 - 31*v(i)*v(i-1)/3 + 
                  4*v(i-2)**2/3 - 19*v(i-2)*v(i-1)/3 + 25*v(i-1)**2/3)

def solve_sos_decomposition(expanded_expr, variables):
    """
    求解形如 c1*L1^2 + c2*L2^2 的SOS分解
    
    参数:
        expanded_expr: 展开后的表达式
        variables: 变量列表，如 [v(i), v(i+1), v(i+2)]
    
    返回:
        分解结果列表 [(c1, L1), (c2, L2)]
    """
    print(f"\n{'='*80}")
    print(f"求解变量的SOS分解: {variables}")
    print(f"{'='*80}")
    
    # 1. 构造对称矩阵Q，使得 expr = v^T * Q * v
    n = len(variables)
    Q = sp.zeros(n, n)
    
    print("\n步骤1: 构造二次型矩阵 Q")
    print("-" * 40)
    
    # 对角线元素
    for i in range(n):
        coeff = expanded_expr.coeff(variables[i]**2)
        Q[i, i] = coeff
        print(f"Q[{i},{i}] = {variables[i]}² 的系数 = {coeff}")
    
    # 非对角线元素（表达式中 2*x*y 的系数对应 Q[i,j] + Q[j,i] = 2*Q[i,j]）
    for i in range(n):
        for j in range(i+1, n):
            coeff = expanded_expr.coeff(variables[i]*variables[j])
            # 由于表达式中是 var_i*var_j，矩阵中对应 2*Q[i,j]
            Q[i, j] = coeff / 2
            Q[j, i] = coeff / 2
            print(f"Q[{i},{j}] = Q[{j},{i}] = {variables[i]}*{variables[j]} 系数/2 = {coeff}/2 = {coeff/2}")
    
    print(f"\n得到的对称矩阵 Q:")
    sp.pprint(Q)
    
    # 2. 特征值分解（理论上有两个非零特征值）
    print("\n步骤2: 矩阵的特征值分解")
    print("-" * 40)
    eigenvals = Q.eigenvals()
    print("特征值及其重数:")
    for val, mult in eigenvals.items():
        print(f"  λ = {sp.nsimplify(val)} (重数: {mult})")
    
    # 3. 寻找秩-1分解
    # 理论上 Q = c1*w1*w1^T + c2*w2*w2^T
    print("\n步骤3: 求解秩-1分解")
    print("-" * 40)
    
    # 方法：通过比较法直接求解
    # 设未知线性表达式: L1 = a1*var0 + a2*var1 + a3*var2
    # 设未知线性表达式: L2 = b1*var0 + b2*var1 + b3*var2
    # 则 Q = c1*[a1,a2,a3]^T*[a1,a2,a3] + c2*[b1,b2,b3]^T*[b1,b2,b3]
    
    # 未知系数
    c1, c2 = sp.symbols('c1 c2')
    a1, a2, a3 = sp.symbols('a1 a2 a3')
    b1, b2, b3 = sp.symbols('b1 b2 b3')
    
    # 构造秩-1矩阵
    w1 = sp.Matrix([a1, a2, a3])
    w2 = sp.Matrix([b1, b2, b3])
    R1 = c1 * (w1 * w1.T)
    R2 = c2 * (w2 * w2.T)
    R = R1 + R2
    
    print("假设分解形式: Q = c1*[a1,a2,a3]ᵀ[a1,a2,a3] + c2*[b1,b2,b3]ᵀ[b1,b2,b3]")
    
    # 4. 比较矩阵元素，建立方程组
    equations = []
    for i_idx in range(n):
        for j_idx in range(n):
            equations.append(sp.Eq(R[i_idx, j_idx], Q[i_idx, j_idx]))
    
    print(f"\n建立 {len(equations)} 个方程:")
    for idx, eq in enumerate(equations[:6]):  # 显示前6个
        print(f"  方程 {idx+1}: {eq}")
    if len(equations) > 6:
        print(f"  ... 还有 {len(equations)-6} 个方程")
    
    # 5. 添加约束条件求解
    # 约束1: 系数是有理数
    # 约束2: 消除缩放歧义（令某些系数为特定值）
    # 约束3: c1, c2 应为正数（因为是平方和）
    
    print("\n步骤4: 求解方程组（添加约束消除缩放歧义）")
    print("-" * 40)
    print("添加约束: b1=1（固定第一个系数）")
    
    # 实际求解时使用数值方法先找到近似解
    equations_constrained = equations + [sp.Eq(b1, 1)]
    
    # 使用 nsolve 寻找数值解（需要提供初始值）
    print("\n使用数值求解作为初始猜测...")
    try:
        # 初始猜测
        initial_guess = {
            a1: 1, a2: -2, a3: 1,
            b1: 1, b2: -2, b3: 1,
            c1: 1, c2: 1
        }
        
        # 选择9个独立方程求解9个未知数
        sol = sp.nsolve(equations_constrained[:9], 
                       [a1, a2, a3, b1, b2, b3, c1, c2], 
                       [1, -2, 1, 1, -2, 1, 1, 1],
                       tol=1e-14, maxsteps=100)
        
        print(f"数值解: {sol}")
        
        # 根据数值解模式，推断符号解
        # 观察数值解的模式，手动构造符号解
        
    except Exception as e:
        print(f"数值求解遇到错误: {e}")
        print("切换到符号模式匹配...")
    
    # 6. 使用启发式方法：寻找整数/有理数解
    # 观察矩阵Q的结构，尝试匹配模式
    
    print("\n步骤5: 启发式模式匹配")
    print("-" * 40)
    
    # 对于beta0，观察系数模式
    # Q = [[10/3, -31/6, 11/6],
    #      [-31/6, 25/3, -19/6],
    #      [11/6, -19/6, 4/3]]
    
    # 尝试寻找形如 (1, -2, 1) 的模式
    # 这是差分算子的典型模式
    
    test_vector1 = sp.Matrix([1, -2, 1])  # 二阶差分
    test_vector2 = sp.Matrix([1, -1, 0])  # 一阶差分
    test_vector3 = sp.Matrix([3, -4, 1])  # WENO的特定组合
    
    # 测试哪个向量能匹配
    print("测试候选向量:")
    
    candidates = [
        ("[1, -2, 1] (二阶差分)", sp.Matrix([1, -2, 1])),
        ("[1, -1, 0]", sp.Matrix([1, -1, 0])),
        ("[3, -4, 1] (WENO专用)", sp.Matrix([3, -4, 1])),
    ]
    
    for name, vec in candidates:
        # 计算 rank(Q - λ*vv^T)
        vvt = vec * vec.T
        print(f"\n  测试向量 {name}:")
        sp.pprint(vec.T)
        
        # 尝试用最小二乘拟合求解λ
        # 对于矩阵的每个非零元素位置
        ratios = []
        for ii in range(n):
            for jj in range(n):
                if vvt[ii, jj] != 0:
                    ratios.append(Q[ii, jj] / vvt[ii, jj])
        
        if ratios:
            avg_ratio = sum(ratios) / len(ratios)
            print(f"    平均比例系数: {sp.nsimplify(avg_ratio)}")
            
            # 验证剩余矩阵
            residual = Q - sp.nsimplify(avg_ratio) * vvt
            print(f"    剩余矩阵的秩: {residual.rank()}")
            
            if residual.rank() == 1:
                print(f"    ✓ 成功！剩余矩阵秩为1，可以继续分解")
                # 寻找第二个向量
                # 剩余矩阵应为 c2 * w2 * w2^T
                for name2, vec2 in candidates:
                    if name2 != name:
                        vvt2 = vec2 * vec2.T
                        ratios2 = []
                        for ii in range(n):
                            for jj in range(n):
                                if vvt2[ii, jj] != 0 and residual[ii, jj] != 0:
                                    ratios2.append(residual[ii, jj] / vvt2[ii, jj])
                        if ratios2:
                            avg_ratio2 = sum(ratios2) / len(ratios2)
                            residual2 = residual - sp.nsimplify(avg_ratio2) * vvt2
                            if residual2.norm() < 1e-10:
                                print(f"    第二个向量 {name2}, 系数: {sp.nsimplify(avg_ratio2)}")
                                return [(sp.nsimplify(avg_ratio), vec), (sp.nsimplify(avg_ratio2), vec2)]

    # 如果没有找到完美匹配，使用特征值方法
    print("\n  使用特征值分解方法...")
    eigenvecs = Q.eigenvects()
    
    sos_decomp = []
    for val_mult in eigenvecs:
        val = sp.nsimplify(val_mult[0])
        mult = val_mult[1]
        if abs(val) > 1e-10:
            # 获取特征向量
            vec = val_mult[2][0]  # 第一个特征向量
            # 归一化
            vec_simplified = sp.nsimplify(vec / sp.sqrt(vec.dot(vec)))
            sos_decomp.append((val, vec_simplified))
    
    return sos_decomp

def format_sos_result(decomp, variables):
    """格式化SOS分解结果为可读形式"""
    if not decomp:
        return "未找到分解"
    
    result = []
    for idx, (coeff, vec) in enumerate(decomp):
        # 构建线性表达式
        terms = []
        for var_idx, var in enumerate(variables):
            if vec[var_idx] != 0:
                term = sp.nsimplify(vec[var_idx]) * var
                terms.append(str(term))
        
        linear_expr = " + ".join(terms).replace("+ -", "- ")
        result.append(f"  项{idx+1}: {sp.nsimplify(coeff)} × ({linear_expr})²")
    
    return "\n".join(result)

# 对每个beta进行分解
print("\n" + "="*80)
print("SOS分解逆向求解")
print("="*80)

# β0分解
vars0 = [v(i), v(i+1), v(i+2)]
decomp0 = solve_sos_decomposition(beta0_expanded, vars0)
print("\n最终结果:")
print(format_sos_result(decomp0, vars0))

# β1分解
vars1 = [v(i-1), v(i), v(i+1)]
decomp1 = solve_sos_decomposition(beta1_expanded, vars1)
print("\nβ1的SOS分解:")
print(format_sos_result(decomp1, vars1))

# β2分解
vars2 = [v(i-2), v(i-1), v(i)]
decomp2 = solve_sos_decomposition(beta2_expanded, vars2)
print("\nβ2的SOS分解:")
print(format_sos_result(decomp2, vars2))

# 验证分解
print("\n="*80)
print("验证分解正确性")
print("="*80)

def verify_decomposition(original_expanded, decomp, variables):
    """验证分解是否正确重构原式"""
    if not decomp:
        return False
    
    reconstructed = 0
    for coeff, vec in decomp:
        # 构造线性表达式
        linear = sum(sp.nsimplify(vec[idx]) * var for idx, var in enumerate(variables))
        reconstructed += sp.nsimplify(coeff) * linear**2
    
    reconstructed = sp.expand(reconstructed)
    diff = sp.simplify(reconstructed - original_expanded)
    return diff == 0

# 验证每个分解
for idx, (decomp, vars_list, original) in enumerate([
    (decomp0, vars0, beta0_expanded),
    (decomp1, vars1, beta1_expanded),
    (decomp2, vars2, beta2_expanded)
], 0):
    print(f"\nβ{idx} 验证:")
    is_correct = verify_decomposition(original, decomp, vars_list)
    print(f"  分解正确: {is_correct}")
    
    if is_correct:
        print(f"  重建表达式: {sp.expand(sum(sp.nsimplify(coeff)*sum(sp.nsimplify(vec[idx])*var for idx, var in enumerate(vars_list))**2 for coeff, vec in decomp))}")