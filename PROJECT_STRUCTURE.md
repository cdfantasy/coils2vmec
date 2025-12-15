# Coils2VMEC 项目结构

## 目录组织

```
coils2vmec/
├── src/                          # 所有源代码
│   └── coils2vmec/               # Python 包
│       ├── __init__.py          # 包初始化和导出
│       ├── fieldline.py         # 磁场线追踪
│       ├── iota.py              # 旋转变换计算
│       ├── lcfs.py              # LCFS 检测
│       ├── descur.py            # DESCUR 接口
│       ├── descur_python.py     # DESCUR 实现（本地模块）
│       ├── plotting.py          # 可视化函数
│       ├── utils.py             # 工具函数
│       ├── fieldline_tracer.py  # f90wrap 生成的接口
│       ├── _fieldline_tracer.so # 编译的 Fortran 扩展
│       └── fortran/             # Fortran 源代码
│           ├── fieldline_tracer_module.f90
│           ├── DLSODE.f
│           ├── hybrd.f
│           └── traceline.f90
├── examples/                     # 示例代码
│   └── example_pwO.py
├── test/                         # 测试数据
├── deprecated/                   # 已弃用文件
│   ├── coils2vmec.py            # 旧的单文件版本
│   └── descur_python.py         # 根目录的副本
├── setup.py                      # 安装脚本（处理 Fortran 编译）
├── pyproject.toml               # PEP 517 配置
├── makefile                      # Fortran 编译 Makefile
├── requirements.txt              # 运行时依赖
├── README.md                     # 项目说明
└── INSTALL.md                    # 安装指南
```

## 文件说明

### 已删除的文件
- ❌ `setup_final.py` - 不需要，setup.py 已经足够
- ❌ `setup_new.py` - 不需要，使用 setup.py

### 源码位置
- ✅ **Python 包**: `src/coils2vmec/` (所有 .py 文件)
- ✅ **Fortran 源码**: `src/coils2vmec/fortran/` (所有 .f90, .f 文件)
- ✅ **编译产物**: `src/coils2vmec/` (.so 和生成的 .py 文件)

## 安装方法

```bash
# 一键安装（推荐）
pip install -e .

# 手动编译 + 安装
make clean && make build
pip install -e .
```

## 导入方式

```python
# 导入整个包
import coils2vmec as c2v

# 使用导出的函数
axis_rz = c2v.find_axis([1.32, 0.0])
lcfs_rz = c2v.find_lcfs([1.32, 0.0])

# 访问 Fortran 模块（如果需要）
from coils2vmec import fieldline_tracer
```

## 关键特性

1. **模块化设计**: 每个功能独立模块 1. **清晰的目录结构**: src/ 用于源码，fortran/ 用于 Fortran 文件
3. **PEP 517 兼容**: 使用 pyproject.toml 标准配置
4. **自动编译**: pip install 时自动处理 Fortran 编译
5. **可编辑安装**: -e 模式便于开发调试

