# Radial Distribution Function
Radial Distribution Function(rdf) 镜像分布函数计算


## 参考
[[VMD] 径向分布函数求配位数的问题求助@sobereva](http://bbs.keinsci.com/thread-15102-1-1.html)<br>
[9041 周期性 RDF 计算程序 CryRDF](https://zhuanlan.zhihu.com/p/178610319)<br>


## 径向分布函数rdf定义
以A原子为中心,半径为$r$的球内,共计有$N_{AB}(r)$个B原子

$$
N_{AB}(r)= \rho_B \int  g_{AB}(r) 4 \pi r^2  dr
$$

其中$\rho_B$是B原子的密度, 为常数. 当$r \to \infty$应有关系

$$
\lim_{r \to \infty} N_{AB}(r) = \rho_B  \frac{4}{3} \pi r^3 = \rho_B \int   4 \pi r^2  dr
$$

所以

$$
\lim_{r \to \infty}  g_{AB}(r) == 1
$$

可得$g(r)$的定义式

$$
g(r)=\lim_{dr \to 0} \frac{dN_{AB}(r)}{\rho_B 4\pi r^2 dr } \simeq  \frac{\Delta N}{\rho_B \Delta V}
$$

所以$\rho_B g(r)$就是以A为球心, B原子在$r$处的密度


## 程序实现

选中第i个A原子坐标`A[i,0:3]`,计算所有B原子与A原子的距离`dis=|B[:,0:3]-A[i,0:3]|`<br>
计算在
$$
\Delta V(r) = V(r \to r+dr )=\frac{4}{3}\pi \left ( (r+dr)^3 -r^3 \right ) = \frac{4}{3} \pi dr \left ( 3r^2 + 3rdr +dr^2 \right )  
$$

内的B原子数$\Delta N(r)= N(r \to r+dr)$, 即`dis[ (dis >= r) and (dis < r+dr) ]`的元素数量.<br>
在程序中,可以使用整除向下取余的方式计算.<br>
- 间距$dr$= `dr`, 区间`[minr,maxr]`,取点数`ngrid=int((maxr-minr)/dr)`
- r的采样网格点$r_{j=0,ngrid}$=`grid=[0,1,2,3,...,ngrid]*dr+minr`
- 计算$B_j$和$A_i$原子距离`dis[j]=|B[j,0:3]-A[i,0:3]|`所处区间`j=floor(dis[j]-minr)/dr)`,<br>
   这里采用floor向下取整表示`(dis >= r) and (dis < r+dr)`<br>
   $N(r_j \to r_j+dr) = N(r_j \to r_j+dr) + 1$
- 统计完所有距离$A_i$小于`maxr`的B原子后,计算$g(r_j)=\frac{N(r_j \to r_j+dr)}{\rho_B \Delta V(r_j)}$
- 可以外套一层循环以所有的A原子为中心计算一遍,最后处以A原子数平均
- **注意,在成键处的$r_b$,因为键长固定,因此$g(r_B)$随着dr反比变化**

## 代码
仓库地址 [Radial_Distribution_Function_rdf](https://github.com/cndaqiang/Radial_Distribution_Function_rdf)<br>
Start-fork-download<br>
![](/uploads/2020/11/rdf.png)

环境
- python3
- Module: numpy, matplotlib

使用示例
- 查看帮助
```
(python37) cndaqiang@mac Desktop$ ./rdf.py
Usage: ./rdf.py [ xxx.xyz | xxx.vasp ] [ none | a | a b c]
```
- 计算晶格常数为30Ang的立方晶胞的rdf
```
(python37) cndaqiang@mac Desktop$ ./rdf.py 303030.xyz 30
read from 303030.xyz
XYZ: set cell to
a: [30.  0.  0.]
b: [ 0. 30.  0.]
c: [ 0.  0. 30.]
nstep,ntyp,nat 1 ['O' 'H'] [ 884 1768]
Super cell [1. 1. 1.]
['O-O' 'O-H' 'H-H']
Save rdf to 303030.xyz.O-O.averagerdf.dat
Save rdf to 303030.xyz.O-H.averagerdf.dat
Save rdf to 303030.xyz.H-H.averagerdf.dat
Save png to 303030.xyz.averagerdf.png
```
- 计算晶格常数`a=3.967 b=8.121 c=8.320`的长方型原胞rdf
```
(python37) cndaqiang@mac Desktop$ ./rdf.py 123.xyz 3.967 8.121 8.320
read from 123.xyz
XYZ: set cell to
a: [3.967 0.    0.   ]
b: [0.    8.121 0.   ]
c: [0.   0.   8.32]
nstep,ntyp,nat 1 ['O' 'H'] [12 24]
Super cell [3. 1. 1.]
['O-O' 'O-H' 'H-H']
Save rdf to 123.xyz.O-O.averagerdf.dat
Save rdf to 123.xyz.O-H.averagerdf.dat
Save rdf to 123.xyz.H-H.averagerdf.dat
Save png to 123.xyz.averagerdf.png
```
- vasp格式
```
(python37) cndaqiang@mac Desktop$ ./rdf.py 123.vasp
read from 123.vasp
nstep,ntyp,nat 1 ['H' 'O'] [24 12]
Super cell [3. 1. 1.]
['H-H' 'H-O' 'O-O']
Save rdf to 123.vasp.H-H.averagerdf.dat
Save rdf to 123.vasp.H-O.averagerdf.dat
Save rdf to 123.vasp.O-O.averagerdf.dat
Save png to 123.vasp.averagerdf.png
```

## 结果展示
```
(python37) cndaqiang@mac Desktop$ ./rdf.py cubic.xyz 12 12 12
read from cubic.xyz
XYZ: set cell to
a: [12.  0.  0.]
b: [ 0. 12.  0.]
c: [ 0.  0. 12.]
nstep,ntyp,nat 1 ['O' 'H'] [ 61 122]
Super cell [2. 2. 2.]
['O-O' 'O-H' 'H-H']
Save rdf to cubic.xyz.O-O.averagerdf.dat
Save rdf to cubic.xyz.O-H.averagerdf.dat
Save rdf to cubic.xyz.H-H.averagerdf.dat
Save png to cubic.xyz.averagerdf.png
```
与VMD计算结果对比
![](/uploads/2020/11/vmd.png)
![](/uploads/2020/11/cubic.xyz.averagerdf.png)
