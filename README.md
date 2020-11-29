# Radial Distribution Function
Radial Distribution Function(rdf) 镜像分布函数计算


## 参考
[[VMD] 径向分布函数求配位数的问题求助@sobereva](http://bbs.keinsci.com/thread-15102-1-1.html)<br>
[9041 周期性 RDF 计算程序 CryRDF](https://zhuanlan.zhihu.com/p/178610319)<br>

## 程序内容解释
[镜像分布函数rdf计算方法和程序](https://cndaqiang.github.io//2020/11/29/rdf/)

## 环境
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
![](https://cndaqiang.github.io/uploads/2020/11/vmd.png)
![](https://cndaqiang.github.io/uploads/2020/11/cubic.xyz.averagerdf.png)
