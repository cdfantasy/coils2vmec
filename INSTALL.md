# 安装指南 (Installation Guide)

## 前置要求 (Prerequisites)

### 系统依赖
- **编译器**: gfortran (Fortran 编译器)
- **Python**: 3.8 或更高版本
- **包管理器依赖**: f90wrap, numpy

### Ubuntu/Debian
```bash
sudo apt-get update
sudo apt-get install -y gfortran python3-dev python3-pip
```

### macOS (使用 Homebrew)
```bash
brew install gcc
pip install python3
```

### CentOS/RHEL
```bash
sudo yum install -y gcc-gfortran python3-devel python3-pip
```

## 完整安装流程

### 1. 克隆仓库
```bash
git clone <repository-url>
cd coils2vmec
```

### 2. 创建虚拟环境（推荐）
```bash
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

### 3. 安装构建依赖
```bash
pip install --upgrade pip setuptools wheel
pip install f90wrap numpy scipy
```

### 4. 编译 Fortran 模块

#### 方式 A: 使用 Makefile（推荐，已配置好）
```bash
# 清理之前的构建（如果有）
make clean

# 编译 Fortran 模块
make

# 验证编译成功
python -c "from fieldline_tracer import fieldline_tracer; print('✓ Fortran module compiled successfully')"
```

#### 方式 B: 使用 setup.py（自动编译）
```bash
# 这会自动调用 Makefile 编译 Fortran
pip install -e .
```

### 5. 安装 Python 依赖
```bash
pip install -r requirements.txt
```

### 6. 安装 descur_python（如果还未安装）
```bash
# 如果 descur_python 在 PyPI 上可用
pip install descur_python

# 或者从源代码安装
# git clone <descur_python-repo>
# cd descur_python
# pip install -e .
```

### 7. 验证安装
```bash
python -c "
import coils2vmec as c2v
print('✓ coils2vmec imported successfully')
print('Available functions:', len(c2v.__all__))
from fieldline_tracer import fieldline_tracer
print('✓ Fortran modules ready')
"
```

## 快速开始

### 运行示例
```bash
cd examples
python example_pwO.py
```

### 使用 Python API
```python
import coils2vmec as c2v
import numpy as np

# 读取线圈文件
coils_data = c2v.read_coils_file('path/to/coils.file')

# 追踪磁力线
result = c2v.trace_fieldlines_parallel(
    initial_guess=np.array([1.32, 0.0]),
    n_fieldlines=99,
    nturn=400,
    nphi=360,
    coils_data=coils_data
)

# 保存结果
c2v.save_fieldlines_hdf5(result, 'fieldlines.h5')
```

## Makefile 说明

Makefile 处理 Fortran 编译，包含以下目标：

```bash
make                  # 编译所有 Fortran 模块
make clean           # 清理编译产物
make distclean        # 彻底清理（包括 f90wrap 生成的文件）
```

### Makefile 变量（如需调整）
- `F90WRAP_FLAGS`: f90wrap 标志
- `FFLAGS`: Fortran 编译标志
- `F2PYFLAGS`: f2py 标志

## 故障排除

### 问题 1: "f90wrap not found"
**解决**: 安装 f90wrap
```bash
pip install f90wrap
```

### 问题 2: "gfortran not found"
**解决**: 安装 Fortran 编译器
```bash
# Ubuntu/Debian
sudo apt-get install gfortran

# macOS
brew install gcc

# CentOS/RHEL
sudo yum install gcc-gfortran
```

### 问题 3: "Cannot find module (...)"
**解决**: 确保 Fortran 编译成功
```bash
make clean
make
```

### 问题 4: 导入 fieldline_tracer 时出错
**解决**: 检查 Python 路径
```python
import sys
print(sys.path)
# 确保包含 fieldline_tracer 编译的目录
```

### 问题 5: descur_python 导入失败
**解决**: 检查是否正确安装
```bash
pip list | grep descur
# 如果没有，安装它
pip install descur_python
```

## 更新依赖

如果需要更新到最新的依赖版本：

```bash
pip install --upgrade -r requirements.txt
```

或者针对特定包：
```bash
pip install --upgrade numpy scipy simsopt
```

## 开发模式安装

如果要进行开发并希望代码改动立即生效：

```bash
# 编辑模式安装
pip install -e .

# 现在修改 coils2vmec 中的 Python 代码会立即生效
# 修改 Fortran 代码需要重新运行 make 编译
```

## 卸载

```bash
# 卸载包
pip uninstall coils2vmec

# 清理编译产物
make distclean
```

## 系统集成

### 添加到 PYTHONPATH（可选）
如果不想用虚拟环境，可以添加到系统 Python 路径：

```bash
export PYTHONPATH="${PYTHONPATH}:/path/to/coils2vmec"
```

### 在 ~/.bashrc 或 ~/.zshrc 中持久化
```bash
echo 'export PYTHONPATH="${PYTHONPATH}:/path/to/coils2vmec"' >> ~/.bashrc
source ~/.bashrc
```

## 许可证

[您的许可证]

## 联系方式

如有问题，请提交 Issue 或联系维护者。
