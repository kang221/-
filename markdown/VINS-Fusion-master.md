## VINS-Fusion-master

#### 文件指南针

camera_models：各种相机模型适配

global_fusion：全局优化模块

loop_fusion：回环检测模块

##### vins_estimater：核心状态估计器

--cmake

​	--FindEigen.cmake

CMake的查找模块，能够让CMake自动找到Eigen的include路径，并定义变量供工程使用。



--launch（在ROS项目中，.launch是启动文件，用于一次性启动一个或者多个节点，并设置节点参数、话题和配置。是XML格式）

​	--vins_rviz.launch

用于启动VINS-Fusion算法节点并在RViz中可视化轨迹和地图。

*什么是RViz？*

RViz是ROS里面的一个可视化工具，用来在三维空间里面显示传感器数据、机器人状态、轨迹、地图等信息。本质上是一个三维可视化“窗口”，帮助你直观观察机器人运行和算法结果。

作用：

- 显示传感器数据：多数据源显示，如点云、IMU等。
- 显示机器人状态：位置、姿态等。
- 显示算法结果：SLM/VIO轨迹、地图、障碍物、路径规划等。
- 调试和演示：可以交互式打开/关闭不同显示模块，方便开发。



**--src**

​	--estimator：状态估计器

​	--factor：因子模块

​	--featureTracker：特征跟踪模块

​	--initial：初始化模块

​	--utility：工具函数

--CMakeLists.txt

--package.xml

config：配置文件目录

docker：环境创建目录

support_files：辅助文件