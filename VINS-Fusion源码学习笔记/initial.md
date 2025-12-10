#### initial

共有initial_aligment、initial_ex_rotation、initial_sfm、solve_5pts四个部分。

##### initial_aligment.cpp

该代码用于对齐IMU和视觉观测数据，共包含了5个函数，其中两个是在.h文件中定义的。

solveGyroscopeBias(.h)、TangentBasis、RefineGravity、LinearAlignment、VisualIMUAlignment(.h）

###### solveGyroscopeBias

该函数用于求解陀螺仪的偏差，对应文章陀螺仪偏差校准部分。

代码中涉及到的旋转，计算的是相邻两帧之间的相对旋转（IMU坐标系下），不是相机到IMU坐标系的转化。
$$
\Delta q_{\text{visual}}
= (q_{bi}^v)^{-1} \otimes q_{bj}^v
$$

~~~cpp
void solveGyroscopeBias(map<double, ImageFrame> &all_image_frame, Vector3d* Bgs)
{
    Matrix3d A;
    Vector3d b;
    Vector3d delta_bg;
    A.setZero();
    b.setZero();
    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    
    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++)
    {
        frame_j = next(frame_i);
        MatrixXd tmp_A(3, 3);//局部雅可比矩阵
        tmp_A.setZero();
        VectorXd tmp_b(3);//误差向量
        tmp_b.setZero();
        Eigen::Quaterniond q_ij(frame_i->second.R.transpose() * frame_j->second.R);
        //四元数q_ij,计算两帧旋转的相对旋转矩阵
        tmp_A = frame_j->second.pre_integration->jacobian.template block<3, 3>(O_R, O_BG);
        //取frame_j IMU预积分雅可比矩阵中关于陀螺仪偏置的3*3子块。
        tmp_b = 2 * (frame_j->second.pre_integration->delta_q.inverse() * q_ij).vec();
        //tmp_b是两帧之间的旋转误差，用于最小二乘球陀螺仪偏置
        //根据IMU预计分得到的旋转增量计算 测量旋转和真实旋转直接的误差，然后将四元数误差转化为3*1向量。乘2倍是未来线性化公式。
        A += tmp_A.transpose() * tmp_A;//所有连续帧雅可比累加
        b += tmp_A.transpose() * tmp_b;//所有连续帧误差累加
    }
    
    //最小二乘法求陀螺仪偏置
    //ldlt()是一种高效的矩阵分解方法，用于解对称正定矩阵
    delta_bg = A.ldlt().solve(b);
    ROS_WARN_STREAM("gyroscope bias initial calibration " << delta_bg.transpose());//打印除陀螺仪偏置值

    //更新滑动窗口内的所有陀螺仪的偏置
    for (int i = 0; i <= WINDOW_SIZE; i++)
        Bgs[i] += delta_bg;

    //使用新偏置重新进行IMU预积分
    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end( ); frame_i++)
    {
        frame_j = next(frame_i);
        //next是指向第i帧的下一帧的意思
        frame_j->second.pre_integration->repropagate(Vector3d::Zero(), Bgs[0]);
        //frame_j是一个map<double,ImageFrame>,frame_j->second就是ImageFrame
       //拆开成三个部分来看：1.取map，frame_j->second 2.取对象成员，.pre_intergation 3.调用成员函数->repropagate(Vector3d::Zero(),Bgs[0])
    }
}
~~~



###### TangentBasis

求重力的切线空间基。

```cpp
MatrixXd TangentBasis(Vector3d &g0)
{
    Vector3d b, c;
    Vector3d a = g0.normalized();//归一化重力向量，得到单位重力向量。
    Vector3d tmp(0, 0, 1);//构建一个临时向量。
    if(a == tmp)
        tmp << 1, 0, 0;
    b = (tmp - a * (a.transpose() * tmp)).normalized();
    //计算tmp在方向a上的投影 a * (a^T * tmp)
    //tmp减去这个投影得到了与a垂直的向量，然后归一化
    c = a.cross(b);//叉乘得到第二个正交单位向量c
    MatrixXd bc(3, 2);
    bc.block<3, 1>(0, 0) = b;
    bc.block<3, 1>(0, 1) = c;
    return bc;
}
```



###### RefineGravity

重力向量细化。求出重力的切线空间基之后，带入进行迭代求解，得到更高精度的重力。

在这段代码中，需要注意的是对于IMU的残差，权重矩阵采用的不是真实矩阵而是单位权矩阵。

```cpp
void RefineGravity(map<double, ImageFrame> &all_image_frame, Vector3d &g, VectorXd &x)
{
    //重力初始化
    Vector3d g0 = g.normalized() * G.norm();
    Vector3d lx, ly;
    //VectorXd x;
    int all_frame_count = all_image_frame.size();
    //计算状态向量的维度，每帧优化三个变量+重力方向在切空间的2维自由度+尺度（s）
    int n_state = all_frame_count * 3 + 2 + 1;

    MatrixXd A{n_state, n_state};
    A.setZero();
    VectorXd b{n_state};
    b.setZero();

    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    
    //四次高斯牛顿迭代（经验循环次数，四次足够精确）
    for(int k = 0; k < 4; k++)
    {
        MatrixXd lxly(3, 2);
        lxly = TangentBasis(g0);
        int i = 0;
        for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++, i++)
        {
            frame_j = next(frame_i);

            MatrixXd tmp_A(6, 9);//建立tmp_A矩阵，共六行残差，对9个状态量求偏导。分别是i时刻的速度（3维）、j时刻的速度、重力（2自由度）、尺度因子（s）
            tmp_A.setZero();
            VectorXd tmp_b(6);//六行残差，3位置残差，3速度残差
            tmp_b.setZero();

            double dt = frame_j->second.pre_integration->sum_dt;
            //从frame_j中层层提取出来两帧（i、j）之间经过的时间
            
            //block的用法
            //block<行数,列数>(起始行，起始列)
            
            //下面四行均是作用于位置残差
            tmp_A.block<3, 3>(0, 0) = -dt * Matrix3d::Identity();//位置
            //Matrix3d::Identity() 3x3单位矩阵
            //填充-dt*I
            tmp_A.block<3, 2>(0, 6) = frame_i->second.R.transpose() * dt * dt / 2 * Matrix3d::Identity() * lxly;//加速度计偏置
            //Ri^T*dt*dt*I*lxly/2；；
            tmp_A.block<3, 1>(0, 8) = frame_i->second.R.transpose() * (frame_j->second.T - frame_i->second.T) / 100.0;//尺度因子
            //Ri^T*(Tj-Ti)/100.0
            tmp_b.block<3, 1>(0, 0) = frame_j->second.pre_integration->delta_p + frame_i->second.R.transpose() * frame_j->second.R * TIC[0] - TIC[0] - frame_i->second.R.transpose() * dt * dt / 2 * g0;
            //位置残差：delta_p+Ri^T*Rj*TIC[0]-TIC[0]-Ri^T*dt^2*g0/2

            //作用于速度残差部分
            tmp_A.block<3, 3>(3, 0) = -Matrix3d::Identity();
            tmp_A.block<3, 3>(3, 3) = frame_i->second.R.transpose() * frame_j->second.R;
            tmp_A.block<3, 2>(3, 6) = frame_i->second.R.transpose() * dt * Matrix3d::Identity() * lxly;
            tmp_b.block<3, 1>(3, 0) = frame_j->second.pre_integration->delta_v - frame_i->second.R.transpose() * dt * Matrix3d::Identity() * g0;

            //此处采用的是单位矩阵作为权重矩阵来固定IMU残差，虽然不符理论，但是在实际工程中是有效的。
            Matrix<double, 6, 6> cov_inv = Matrix<double, 6, 6>::Zero();
            //cov.block<6, 6>(0, 0) = IMU_cov[i + 1];
            //MatrixXd cov_inv = cov.inverse();
            cov_inv.setIdentity();

            MatrixXd r_A = tmp_A.transpose() * cov_inv * tmp_A;
            VectorXd r_b = tmp_A.transpose() * cov_inv * tmp_b;

            //对全局矩阵H进行参数的填充和累加
            //当前帧的IMU约束如何影响自身的位姿和速度等局部参数
            A.block<6, 6>(i * 3, i * 3) += r_A.topLeftCorner<6, 6>();
            b.segment<6>(i * 3) += r_b.head<6>();

            //当前帧的IMU约束如何共同影响尺度或者重力等全局参数
            A.bottomRightCorner<3, 3>() += r_A.bottomRightCorner<3, 3>();
            b.tail<3>() += r_b.tail<3>();

            //说明调整全局状态量会如何影响局部的状态
            A.block<6, 3>(i * 3, n_state - 3) += r_A.topRightCorner<6, 3>();
            A.block<3, 6>(n_state - 3, i * 3) += r_A.bottomLeftCorner<3, 6>();
        }
            A = A * 1000.0;
            b = b * 1000.0;
            x = A.ldlt().solve(b);
            VectorXd dg = x.segment<2>(n_state - 3);
            g0 = (g0 + lxly * dg).normalized() * G.norm();
            //double s = x(n_state - 1);
    }   
    g = g0;
}
```



###### LinearAlignment

求解尺度、重力向量和速度。线性求解得到速度、尺度、重力的初始值，然后再进行重力细化处理，得到更精确的重力值

```cpp
bool LinearAlignment(map<double, ImageFrame> &all_image_frame, Vector3d &g, VectorXd &x)
{
    int all_frame_count = all_image_frame.size();
    //三维速度+三维重力+尺度 
    int n_state = all_frame_count * 3 + 3 + 1;

    MatrixXd A{n_state, n_state};
    A.setZero();
    VectorXd b{n_state};
    b.setZero();

    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    int i = 0;
    
    //遍历每一对关键帧
    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++, i++)
    {
        frame_j = next(frame_i);

        MatrixXd tmp_A(6, 10);
        tmp_A.setZero();
        VectorXd tmp_b(6);
        tmp_b.setZero();

        double dt = frame_j->second.pre_integration->sum_dt;//求相邻关键帧的时间间隔

        tmp_A.block<3, 3>(0, 0) = -dt * Matrix3d::Identity();
        tmp_A.block<3, 3>(0, 6) = frame_i->second.R.transpose() * dt * dt / 2 * Matrix3d::Identity();
        tmp_A.block<3, 1>(0, 9) = frame_i->second.R.transpose() * (frame_j->second.T - frame_i->second.T) / 100.0;     
        tmp_b.block<3, 1>(0, 0) = frame_j->second.pre_integration->delta_p + frame_i->second.R.transpose() * frame_j->second.R * TIC[0] - TIC[0];
        //cout << "delta_p   " << frame_j->second.pre_integration->delta_p.transpose() << endl;
        tmp_A.block<3, 3>(3, 0) = -Matrix3d::Identity();
        tmp_A.block<3, 3>(3, 3) = frame_i->second.R.transpose() * frame_j->second.R;
        tmp_A.block<3, 3>(3, 6) = frame_i->second.R.transpose() * dt * Matrix3d::Identity();
        tmp_b.block<3, 1>(3, 0) = frame_j->second.pre_integration->delta_v;
        //cout << "delta_v   " << frame_j->second.pre_integration->delta_v.transpose() << endl;

        Matrix<double, 6, 6> cov_inv = Matrix<double, 6, 6>::Zero();
        //cov.block<6, 6>(0, 0) = IMU_cov[i + 1];
        //MatrixXd cov_inv = cov.inverse();
        cov_inv.setIdentity();

        MatrixXd r_A = tmp_A.transpose() * cov_inv * tmp_A;
        VectorXd r_b = tmp_A.transpose() * cov_inv * tmp_b;

        A.block<6, 6>(i * 3, i * 3) += r_A.topLeftCorner<6, 6>();
        b.segment<6>(i * 3) += r_b.head<6>();

        A.bottomRightCorner<4, 4>() += r_A.bottomRightCorner<4, 4>();
        b.tail<4>() += r_b.tail<4>();

        A.block<6, 4>(i * 3, n_state - 4) += r_A.topRightCorner<6, 4>();
        A.block<4, 6>(n_state - 4, i * 3) += r_A.bottomLeftCorner<4, 6>();
    }
    A = A * 1000.0;
    b = b * 1000.0;
    x = A.ldlt().solve(b);//求线性方程组
    double s = x(n_state - 1) / 100.0;
    ROS_DEBUG("estimated scale: %f", s);
    g = x.segment<3>(n_state - 4);
    ROS_DEBUG_STREAM(" result g     " << g.norm() << " " << g.transpose());
    if(fabs(g.norm() - G.norm()) > 0.5 || s < 0)
    {
        return false;
    }

    RefineGravity(all_image_frame, g, x);
    s = (x.tail<1>())(0) / 100.0;
    (x.tail<1>())(0) = s;
    ROS_DEBUG_STREAM(" refine     " << g.norm() << " " << g.transpose());
    if(s < 0.0 )
        return false;   
    else
        return true;
}
```



###### VisualIMUAlignment

速度、重力、尺度初始化函数

```cpp
bool VisualIMUAlignment(map<double, ImageFrame> &all_image_frame, Vector3d* Bgs, Vector3d &g, VectorXd &x)
{
    solveGyroscopeBias(all_image_frame, Bgs);

    if(LinearAlignment(all_image_frame, g, x))
        return true;
    else 
        return false;
}
```



##### initial_aligment.h

~~~cpp
#pragma once
//防止头文件被重复包含
#include <eigen3/Eigen/Dense>
//矩阵库
#include <iostream>
//输入输出流
#include "../factor/imu_factor.h"
//IMU因子、用于非线性优化
#include "../utility/utility.h"
//工具函数
#include <ros/ros.h>
//ROS基础功能
#include <map>
#include "../estimator/feature_manager.h"
//图像特征点管理器
~~~

ImageFrame类

​	该类定义了每帧图像的数据结构，包含了时间戳，特征点、位姿、IMU预计分对象还有是否为关键帧。

~~~cpp
class ImageFrame
{
    pubilc:
    	ImageFrame(){};//默认构造函数,无参。
        ImageFrame(const map<int,vector<pair<int, Eigen::Martix<double, 7, 1>>>>& _points, double _t):t{_t},is_key_frame{false}
    {
        points = _points;
    };
    /*带参数的构造函数，共两个参数
     参数一：cost map<int,vector<pair<int, Eigen::Matrix<double, 7 ,1>>>>& _points
     _points是一个键对，（帧：ID，值：特征点列表向量）。
     参数二：double _t
     _t是这帧图像的时间戳。
    */
    	map<int, vecotr<pair<int, Eigen::Matrix<double, 7, 1>> > > points;//定义points存储每帧图像的特征点信息
    double t;//定义t表示图像帧的采集时间
    Matrix3d R;//定义R作为当前帧相机在世界坐标系下的旋转矩阵
    Vector3d T;//定义T作为当前帧在世界坐标系下的位置
    IntegrationBase *pre_integration;//定义指针pre_integraion为IMU预积分对象
    bool is_key_frame;//使用布尔变量，标记当前帧是否是关键帧
}
~~~



##### initial_ex_rotation.cpp

共有CalibrationExRotation（h)、solveRelativeR（h）、testTriangulation（h）、decomposeE（h）

###### InitialEXRotation

类中构造函数的实现

```cpp
InitialEXRotation::InitialEXRotation(){
    frame_count = 0;
    Rc.push_back(Matrix3d::Identity());//相对相机旋转,两关键帧在相机坐标系下的相对旋转。通过sfm得到。
    Rc_g.push_back(Matrix3d::Identity());//通过IMU测量和外参估计间接计算的相机相对旋转。
    Rimu.push_back(Matrix3d::Identity());//IMU相对旋转，两关键帧在IMU坐标系下的相对旋转
    ric = Matrix3d::Identity();//相机到IMU的外参旋转矩阵
}
//Matrix3d是Eigen库中的一个特定类型，表3*3的double矩阵
```



###### CalibrationExRotation

**四元数乘法**转化为线性代数形式，本质是利用四元数的特殊结构，将四元数的复杂运算，转化为4*4的矩阵乘法操作。	

四元数q可以表示为：
$$
\mathbf{q} = w + x\mathbf{i} + y\mathbf{j} + z\mathbf{k}
$$
但是在VINS的应用中，包括Eigen库等，都是:
$$
\mathbf{q} = \begin{pmatrix} \mathbf{v} \\ w \end{pmatrix}
$$
两个四元数 $\mathbf{q}_a = (w_a, \mathbf{v}_a)$和$\mathbf{q}_b = (w_b, \mathbf{v}_b)$ 的乘法$\mathbf{q}_a \otimes \mathbf{q}_b$ 定义为：

$$
\mathbf{q}_a \otimes \mathbf{q}_b = (w_a w_b - \mathbf{v}_a \cdot \mathbf{v}_b) + (w_a \mathbf{v}_b + w_b \mathbf{v}_a + \mathbf{v}_a \times \mathbf{v}_b)
$$

这里 $\cdot$ 是点积，$\times$ 是叉积。

**线性化**：任何四元数乘法$\mathbf{q}_a \otimes \mathbf{q}_b = \mathbf{q}_c$ 都可以通过矩阵乘法表示为
$$
\mathbf{q}_a \otimes \mathbf{q}_b = \mathbf{L}(\mathbf{q}_a) \mathbf{q}_b = \mathbf{R}(\mathbf{q}_b) \mathbf{q}_a
$$
点积的线性化：$\mathbf{v}_a \cdot \mathbf{v}_b = \mathbf{v}_a^T \mathbf{v}_b$

叉积的线性化：$\mathbf{v}_a \times \mathbf{v}_b = [\mathbf{v}_a]_{\times} \mathbf{v}_b$ 其中 $[\mathbf{v}_a]_{\times}$ 是向量 $\mathbf{v}_a$ 的**反对称矩阵**。

$$
\left[\mathbf{v}_a\right]_{\times} = \begin{pmatrix} 0 & -v_{a3} & v_{a2} \\ v_{a3} & 0 & -v_{a1} \\ -v_{a2} & v_{a1} & 0 \end{pmatrix}
$$
根据点积和叉积构造出来的左乘矩阵如下：
$$
\mathbf{L}(\mathbf{q}_a) = \begin{pmatrix} w_a \mathbf{I} + [\mathbf{v}_a]_{\times} & \mathbf{v}_a \\ -\mathbf{v}_a^\top & w_a \end{pmatrix}
$$
同理，右乘矩阵如下：
$$
\mathbf{R}(\mathbf{q}_b) = \begin{pmatrix} w_b \mathbf{I} - [\mathbf{v}_b]_{\times} & \mathbf{v}_b \\ -\mathbf{v}_b^\top & w_b \end{pmatrix}
$$
用于标定相机于IMU之间的旋转外参。


```cpp
bool InitialEXRotation::CalibrationExRotation(vector<pair<Vector3d, Vector3d>> corres, Quaterniond delta_q_imu, Matrix3d &calib_ric_result)
{
    frame_count++;//累加关键帧数量
    Rc.push_back(solveRelativeR(corres));//求两关键帧之间的相对旋转（相机坐标系下）
    Rimu.push_back(delta_q_imu.toRotationMatrix());
    //toRotationMatrix是四元数类的方法，能够将四元数形式，转化为矩阵形式
    Rc_g.push_back(ric.inverse() * delta_q_imu * ric);
    //将IMU测量的相对旋转矩阵通过外参旋转转化到相机坐标系下

    Eigen::MatrixXd A(frame_count * 4, 4);
    //MatrixXd X表示矩阵的行数是动态的，d代表双精度浮点数
    A.setZero();
    int sum_ok = 0;
    for (int i = 1; i <= frame_count; i++)
    {
        Quaterniond r1(Rc[i]);
        Quaterniond r2(Rc_g[i]);

        double angular_distance = 180 / M_PI * r1.angularDistance(r2);
        //angularDistance用于计算两个四元数自己的最短距离
        ROS_DEBUG(
            "%d %f", i, angular_distance);
        //ROS系统中的输出

        //huber是Huber损失函数所派生出的权重因子
        double huber = angular_distance > 5.0 ? 5.0 / angular_distance : 1.0;
        ++sum_ok;
        Matrix4d L, R;
        //L矩阵是左乘矩阵

        double w = Quaterniond(Rc[i]).w();
        //.w 提取四元素标量部分
        //q=(w,q),w是变量，q是向量
        Vector3d q = Quaterniond(Rc[i]).vec();
        //.vec 提取向量部分
        
        //填充L左乘矩阵，目的是为了将四元数乘法转化为线性乘法，简化求解
        L.block<3, 3>(0, 0) = w * Matrix3d::Identity() + Utility::skewSymmetric(q);
        L.block<3, 1>(0, 3) = q;
        L.block<1, 3>(3, 0) = -q.transpose();
        L(3, 3) = w;

        Quaterniond R_ij(Rimu[i]);
        w = R_ij.w();
        q = R_ij.vec();
        R.block<3, 3>(0, 0) = w * Matrix3d::Identity() - Utility::skewSymmetric(q);
        R.block<3, 1>(0, 3) = q;
        R.block<1, 3>(3, 0) = -q.transpose();
        R(3, 3) = w;

        A.block<4, 4>((i - 1) * 4, 0) = huber * (L - R);
    }//将四元数旋转约束乘以鲁棒权重后，堆叠到全局矩阵A的指定位置。

    JacobiSVD<MatrixXd> svd(A, ComputeFullU | ComputeFullV);
    //jacobisvds 是eigen库中实现SVD算法的类。
    //A是第一个参数，进行奇异值求解的矩阵
    //ComputeFullU | ComputeFullV是一个flag,表示要求计算完整的左/右奇异向量矩阵
    Matrix<double, 4, 1> x = svd.matrixV().col(3);
    //提取出SVD右奇异向量矩阵的第四列，第四列是A矩阵的零空间解，也就是最优的相机-IMU外参四元数qic。IMU到相机的选择
    Quaterniond estimated_R(x);
    ric = estimated_R.toRotationMatrix().inverse();//转化为相机到IMU的旋转
    //cout << svd.singularValues().transpose() << endl;
    //cout << ric << endl;
    Vector3d ric_cov;
    ric_cov = svd.singularValues().tail<3>();
    //从奇异值中提取出跟外参Ric估计不确定度相关的三个值。
    
    //满足一下俩条件则认为外参标定成功。
    //条件一：约束数量充足
    //条件二：估计精度达到要求，0.25被认为是一个经验性的阈值检查结果
    if (frame_count >= WINDOW_SIZE && ric_cov(1) > 0.25)
    {
        calib_ric_result = ric;
        return true;
    }
    else
        return false;
}
```



###### solveRelativeR

利用对极几何原理，从视觉匹配的特征点中求解出相邻两帧相机之间的相对旋转矩阵Rc。

OpenCV中，矩阵的存储是行优先，Eigen中是列优先。

```cpp
Matrix3d InitialEXRotation::solveRelativeR(const vector<pair<Vector3d, Vector3d>> &corres)//corres是特征点对的坐标。
{
    if (corres.size() >= 9)
    {
        vector<cv::Point2f> ll, rr;
        //cv是opencv的命名空间，point2f用于表示一个二维浮点坐标点。
        for (int i = 0; i < int(corres.size()); i++)
        {
            ll.push_back(cv::Point2f(corres[i].first(0), corres[i].first(1)));//左图点集合，前一帧
            rr.push_back(cv::Point2f(corres[i].second(0), corres[i].second(1)));//右图点集合，当前帧
        }
        cv::Mat E = cv::findFundamentalMat(ll, rr);
        //根据两相邻关键帧的匹配特征点ll、rr求解基础矩阵F。
        //此处涉及到基础矩阵（描述像素坐标下的对极约束）和本质矩阵（描述归一化坐标下的对极约束）的线代概念，需要学习。
        cv::Mat_<double> R1, R2, t1, t2;
        decomposeE(E, R1, R2, t1, t2);

        //检验求解出的旋转矩阵是否有效。有效的旋转矩阵R必须满足特征值=1
        if (determinant(R1) + 1.0 < 1e-09)
        {
            E = -E;
            decomposeE(E, R1, R2, t1, t2);
        }
        
        //本质矩阵E的分解会产生一下四组解，(R1，t1),(R1,t2),(R2,T1),(R2,t2)
        //但是只有一组能够使得对应的3D点在两个相机中都具有正深度（位于相机前方）
        double ratio1 = max(testTriangulation(ll, rr, R1, t1), testTriangulation(ll, rr, R1, t2));//比较旋转R1搭配t1,t2的效果，消除R1对应平移t歧义。 
        double ratio2 = max(testTriangulation(ll, rr, R2, t1), testTriangulation(ll, rr, R2, t2));//比较旋转R2搭配t1,t2的效果，消除R2对应平移t歧义。 
        //三元运算符，如果ration1>ration2，则取R1，否则，取R2
        cv::Mat_<double> ans_R_cv = ratio1 > ratio2 ? R1 : R2;

        //将以opencv的矩阵存储方式转化为eigen的矩阵存储方式
        Matrix3d ans_R_eigen;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                ans_R_eigen(j, i) = ans_R_cv(i, j);
        return ans_R_eigen;
    }
    return Matrix3d::Identity();
}
```



###### testTriangulation

通过三角化恢复3D点，并检查点是否落在两个相机前方。返回“在前方的点的比例”，用于评估外参是否合理。

```cpp
double InitialEXRotation::testTriangulation(const vector<cv::Point2f> &l,
                                          const vector<cv::Point2f> &r,
                                          cv::Mat_<double> R, cv::Mat_<double> t)
{
    cv::Mat pointcloud;
    //左相机的投影矩阵，设为单位阵意味着左相机作为参考坐标系。
    cv::Matx34f P = cv::Matx34f(1, 0, 0, 0,
                                0, 1, 0, 0,
                                0, 0, 1, 0);
    //右相机的投影矩阵，表示相对于左相机的位姿。
    cv::Matx34f P1 = cv::Matx34f(R(0, 0), R(0, 1), R(0, 2), t(0),
                                 R(1, 0), R(1, 1), R(1, 2), t(1),
                                 R(2, 0), R(2, 1), R(2, 2), t(2));
    cv::triangulatePoints(P, P1, l, r, pointcloud);
    //OpenCV内置库，利用两张图像来对图像中的匹配特征点进行三维重建。
    //输入两个投影矩阵（3*4），两个图像点集；返回一个三维点云（4*N的齐次坐标）。
    int front_count = 0;
    for (int i = 0; i < pointcloud.cols; i++)
    {
        
        double normal_factor = pointcloud.col(i).at<float>(3);
        //从三角化得到的齐次坐标点中取出第四个分量w。

        cv::Mat_<double> p_3d_l = cv::Mat(P) * (pointcloud.col(i) / normal_factor);
        cv::Mat_<double> p_3d_r = cv::Mat(P1) * (pointcloud.col(i) / normal_factor);
        //反齐次化之后，再投影到左右相机坐标系
        
        //判断点是否在两个相机前方
        if (p_3d_l(2) > 0 && p_3d_r(2) > 0)
            front_count++;
    }
    ROS_DEBUG("MotionEstimator: %f", 1.0 * front_count / pointcloud.cols);
    return 1.0 * front_count / pointcloud.cols;
}

```



###### decomposeE

实现了对本质矩阵E的分解，是从两张图像的几何约束中提取相对旋转R和相对平移t的标准数学过程。

```cpp
void InitialEXRotation::decomposeE(cv::Mat E,
                                 cv::Mat_<double> &R1, cv::Mat_<double> &R2,
                                 cv::Mat_<double> &t1, cv::Mat_<double> &t2)
{
    cv::SVD svd(E, cv::SVD::MODIFY_A);
    //对本质矩阵E进行奇异值分解（SVD）
    //MODIFY_A意味着允许修改，节省内存；同时也能强制的将秩设为2，确保本质矩阵是有效的。
    //得到三个矩阵，U左奇异向量矩阵，Σ奇异值矩阵、VT右奇异向量矩阵的转置
    cv::Matx33d W(0, -1, 0,
                  1, 0, 0,
                  0, 0, 1);//旋转90度
    cv::Matx33d Wt(0, 1, 0,
                   -1, 0, 0,
                   0, 0, 1);//旋转-90度
    //构造相对旋转R
    R1 = svd.u * cv::Mat(W) * svd.vt;
    R2 = svd.u * cv::Mat(Wt) * svd.vt;
    //构造相对平移t
    t1 = svd.u.col(2);
    t2 = -svd.u.col(2);
}
```



##### initial_ex_rotation.h

```cpp
#pragma once 

#include <vector>
#include "../estimator/parameters.h"
using namespace std;

#include <opencv2/opencv.hpp>

#include <eigen3/Eigen/Dense>
using namespace Eigen;
#include <ros/console.h>

/* This class help you to calibrate extrinsic rotation between imu and camera when your totally don't konw the extrinsic parameter */
class InitialEXRotation
{
public:
	InitialEXRotation();
    bool CalibrationExRotation(vector<pair<Vector3d, Vector3d>> corres, Quaterniond delta_q_imu, Matrix3d &calib_ric_result);
private:
	Matrix3d solveRelativeR(const vector<pair<Vector3d, Vector3d>> &corres);

    double testTriangulation(const vector<cv::Point2f> &l,
                             const vector<cv::Point2f> &r,
                             cv::Mat_<double> R, cv::Mat_<double> t);
    void decomposeE(cv::Mat E,
                    cv::Mat_<double> &R1, cv::Mat_<double> &R2,
                    cv::Mat_<double> &t1, cv::Mat_<double> &t2);
    int frame_count;

    vector< Matrix3d > Rc;
    vector< Matrix3d > Rimu;
    vector< Matrix3d > Rc_g;
    Matrix3d ric;
};
```



##### initial_sfm.h

```cpp
#pragma once 
#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <cstdlib>
#include <deque>
#include <map>
#include <opencv2/core/eigen.hpp>
#include <opencv2/opencv.hpp>
using namespace Eigen;
using namespace std;

//定义了sfm特征结构体
struct SFMFeature
{
    bool state;
    int id;
    vector<pair<int,Vector2d>> observation;
    double position[3];
    double depth;
};

//定义三维点的重投影误差
struct ReprojectionError3D
{
    //构造函数，用于保存观测点
	ReprojectionError3D(double observed_u, double observed_v)
		:observed_u(observed_u), observed_v(observed_v)
		{}

	template <typename T>//表示T是一个类型参数，可以在编译时由编译器替换为具体类型
    //重载operator()这个在结构体中，让这个结构体具备函数功能
	bool operator()(const T* const camera_R, const T* const camera_T, const T* point, T* residuals) const
	{
		T p[3];
		ceres::QuaternionRotatePoint(camera_R, point, p);
        //用四元素旋转一个三位点
		p[0] += camera_T[0]; p[1] += camera_T[1]; p[2] += camera_T[2];
        //三维点的平移变换
		T xp = p[0] / p[2];
    	T yp = p[1] / p[2];
        //透视投影，将三维点投影到二维归一化图像平面
    	residuals[0] = xp - T(observed_u);
    	residuals[1] = yp - T(observed_v);
        //xp,yp是经过旋转、平移、透视投影后得到的预测归一化图像坐标
        //observed_u、v是实际观测的二维像素
    	return true;
	}

    //静态Crate函数，封装了Ceres的残差函数创建逻辑
    //ChatGPT是下面这样讲的，但是暂时不太能理解
    //Create 就是一个“工厂函数”，把 ReprojectionError3D 对象包装成 Ceres 可以直接用的 CostFunction。
	static ceres::CostFunction* Create(const double observed_x,
	                                   const double observed_y) 
	{
	  return (new ceres::AutoDiffCostFunction<
	          ReprojectionError3D, 2, 4, 3, 3>(
	          	new ReprojectionError3D(observed_x,observed_y)));
	}

	double observed_u;
	double observed_v;
};

//Global Structure-from-Motion全局结构与运动恢复
class GlobalSFM
{
public:
	GlobalSFM();
	bool construct(int frame_num, Quaterniond* q, Vector3d* T, int l,
			  const Matrix3d relative_R, const Vector3d relative_T,
			  vector<SFMFeature> &sfm_f, map<int, Vector3d> &sfm_tracked_points);

private:
	bool solveFrameByPnP(Matrix3d &R_initial, Vector3d &P_initial, int i, vector<SFMFeature> &sfm_f);

	void triangulatePoint(Eigen::Matrix<double, 3, 4> &Pose0, Eigen::Matrix<double, 3, 4> &Pose1,
							Vector2d &point0, Vector2d &point1, Vector3d &point_3d);
	void triangulateTwoFrames(int frame0, Eigen::Matrix<double, 3, 4> &Pose0, 
							  int frame1, Eigen::Matrix<double, 3, 4> &Pose1,
							  vector<SFMFeature> &sfm_f);

	int feature_num;
};
```



##### initial_sfm.cpp

一共有四个函数 triangulatePoint、solveFrameByPnP、triangulateTwoFrames、construct 都在头文件中说明。

###### triangulatePoint

目的是从两帧图像的对应特征点，以及两帧的相机位姿，恢复该点在三维空间中的坐标，这就是三角化（triangulation）。

```cpp
void GlobalSFM::triangulatePoint(Eigen::Matrix<double, 3, 4> &Pose0, Eigen::Matrix<double, 3, 4> &Pose1,
						Vector2d &point0, Vector2d &point1, Vector3d &point_3d)
{
    //创建构造矩阵
	Matrix4d design_matrix = Matrix4d::Zero();
	design_matrix.row(0) = point0[0] * Pose0.row(2) - Pose0.row(0);
	design_matrix.row(1) = point0[1] * Pose0.row(2) - Pose0.row(1);
	design_matrix.row(2) = point1[0] * Pose1.row(2) - Pose1.row(0);
	design_matrix.row(3) = point1[1] * Pose1.row(2) - Pose1.row(1);
	Vector4d triangulated_point;
    
    //用SVD求解X（最小二乘法）
	triangulated_point =
		      design_matrix.jacobiSvd(Eigen::ComputeFullV).matrixV().rightCols<1>();
	
    //非齐次化
    point_3d(0) = triangulated_point(0) / triangulated_point(3);
	point_3d(1) = triangulated_point(1) / triangulated_point(3);
	point_3d(2) = triangulated_point(2) / triangulated_point(3);
}
```



###### solveFrameByPnP

PnP（perspective-n-point）问题。利用已知3D点（来自其他关键帧三角化）+当前帧对应的2D像素，求解当前相机的位姿。

三角化得到的3D坐标是在相对坐标系下的，因此要让当前帧加入统一坐标系，就需要求位姿，使得全局一致。根据SLAM的经验来讲，至少需要15个特征点。

```cpp
bool GlobalSFM::solveFrameByPnP(Matrix3d &R_initial, Vector3d &P_initial, int i,
								vector<SFMFeature> &sfm_f)
{
	vector<cv::Point2f> pts_2_vector;//2D像素点
	vector<cv::Point3f> pts_3_vector;//3D空间点
	for (int j = 0; j < feature_num; j++)
	{
		if (sfm_f[j].state != true)
			continue;//判断特侦点是否完成三角化，完成则继续
		Vector2d point2d;
		for (int k = 0; k < (int)sfm_f[j].observation.size(); k++)
		{
			if (sfm_f[j].observation[k].first == i)
			{
				Vector2d img_pts = sfm_f[j].observation[k].second;
           		cv::Point2f pts_2(img_pts(0), img_pts(1));
				pts_2_vector.push_back(pts_2);
				cv::Point3f pts_3(sfm_f[j].position[0], sfm_f[j].position[1], sfm_f[j].position[2]);
				pts_3_vector.push_back(pts_3);
				break;
			}
		}
	}
    //小于10个点，报错；小于15个提示不稳定
	if (int(pts_2_vector.size()) < 15)
	{
		printf("unstable features tracking, please slowly move you device!\n");
		if (int(pts_2_vector.size()) < 10)
			return false;
	}
	cv::Mat r, rvec, t, D, tmp_r;
	cv::eigen2cv(R_initial, tmp_r);
    //eigen和opencv矩阵格式转换
	cv::Rodrigues(tmp_r, rvec);
    //Rodrigues是旋转矩阵与旋转向量相互转换的函数
    //Rodrigues向量=旋转轴方向 x 旋转角度
	cv::eigen2cv(P_initial, t);
	cv::Mat K = (cv::Mat_<double>(3, 3) << 1, 0, 0, 0, 1, 0, 0, 0, 1);
    //定义相机内参矩阵K
	bool pnp_succ;
	pnp_succ = cv::solvePnP(pts_3_vector, pts_2_vector, K, D, rvec, t, 1);//调用PnP函数求解。
	if(!pnp_succ)
	{
		return false;
	}
	cv::Rodrigues(rvec, r);
	//cout << "r " << endl << r << endl;
	MatrixXd R_pnp;
	cv::cv2eigen(r, R_pnp);
	MatrixXd T_pnp;
	cv::cv2eigen(t, T_pnp);
	R_initial = R_pnp;
	P_initial = T_pnp;
	return true;

}
```



##### solve_5pts.h

```cpp
#pragma once

#include <vector>
using namespace std;

#include <opencv2/opencv.hpp>
//#include <opencv2/core/eigen.hpp>
#include <eigen3/Eigen/Dense>
using namespace Eigen;

#include <ros/console.h>

class MotionEstimator
{
  public:

    bool solveRelativeRT(const vector<pair<Vector3d, Vector3d>> &corres, Matrix3d &R, Vector3d &T);

  //一下两个函数最终在源文件中没有被定义，而是被OpenVC库中的函数所替代
  private:
    double testTriangulation(const vector<cv::Point2f> &l,
                             const vector<cv::Point2f> &r,
                             cv::Mat_<double> R, cv::Mat_<double> t);
    void decomposeE(cv::Mat E,
                    cv::Mat_<double> &R1, cv::Mat_<double> &R2,
                    cv::Mat_<double> &t1, cv::Mat_<double> &t2);
};
```



##### solve_5pts.cpp

5点法，从两帧匹配点中求解两帧之间的相对旋转和相对平移方向。

```cpp
//namspace cv{}的作用是通过命名空间cv来替换OpenCV库中的函数
namespace cv {
    void decomposeEssentialMat( InputArray _E, OutputArray _R1, OutputArray _R2, OutputArray _t )
        //此处没有给出常见的参数类型，InputArray/OutputArray是OpenCV定义的类。
    {
        Mat E = _E.getMat().reshape(1, 3);//E还是3*3的矩阵，reshape（通道，行数）。_E是3*3，reshape不改变数据量，仅改变矩阵的排列方式。
        
        CV_Assert(E.cols == 3 && E.rows == 3);
        //断言语句，检查条件是否成立
        
        Mat D, U, Vt;
        SVD::compute(E, D, U, Vt);

        if (determinant(U) < 0) U *= -1.;
        if (determinant(Vt) < 0) Vt *= -1.;

        //定义旋转辅助矩阵W
        Mat W = (Mat_<double>(3, 3) << 0, 1, 0, -1, 0, 0, 0, 0, 1);
        W.convertTo(W, E.type());//保证W和E的类型一致

        //Essential Matrix分解生成旋转矩阵R和平移向量t
        Mat R1, R2, t;
        R1 = U * W * Vt;
        R2 = U * W.t() * Vt;
        t = U.col(2) * 1.0;//只有方向，没有尺度

        R1.copyTo(_R1);
        R2.copyTo(_R2);
        t.copyTo(_t);
    }

    int recoverPose( InputArray E, InputArray _points1, InputArray _points2, InputArray _cameraMatrix,
                         OutputArray _R, OutputArray _t, InputOutputArray _mask)
    {

        Mat points1, points2, cameraMatrix;
        _points1.getMat().convertTo(points1, CV_64F);
        //先从_points1中取出矩阵，然后再将矩阵转化为64位浮点数
        _points2.getMat().convertTo(points2, CV_64F);
        _cameraMatrix.getMat().convertTo(cameraMatrix, CV_64F);

        int npoints = points1.checkVector(2);
        //检查points2是否可以看作N个2D向量
        CV_Assert( npoints >= 0 && points2.checkVector(2) == npoints &&
                                  points1.type() == points2.type());

        CV_Assert(cameraMatrix.rows == 3 && cameraMatrix.cols == 3 && cameraMatrix.channels() == 1);

        if (points1.channels() > 1)
        {
            points1 = points1.reshape(1, npoints);
            points2 = points2.reshape(1, npoints);
        }

        double fx = cameraMatrix.at<double>(0,0);
        double fy = cameraMatrix.at<double>(1,1);
        double cx = cameraMatrix.at<double>(0,2);
        double cy = cameraMatrix.at<double>(1,2);

        //从像素坐标归一化到相机坐标
        points1.col(0) = (points1.col(0) - cx) / fx;
        points2.col(0) = (points2.col(0) - cx) / fx;
        points1.col(1) = (points1.col(1) - cy) / fy;
        points2.col(1) = (points2.col(1) - cy) / fy;

        points1 = points1.t();//求转置
        points2 = points2.t();

        Mat R1, R2, t;
        decomposeEssentialMat(E, R1, R2, t);//分解Essential Matrix E
        Mat P0 = Mat::eye(3, 4, R1.type());
        Mat P1(3, 4, R1.type()), P2(3, 4, R1.type()), P3(3, 4, R1.type()), P4(3, 4, R1.type());
        //构造四种可能的候选矩阵
        P1(Range::all(), Range(0, 3)) = R1 * 1.0; P1.col(3) = t * 1.0;
        P2(Range::all(), Range(0, 3)) = R2 * 1.0; P2.col(3) = t * 1.0;
        P3(Range::all(), Range(0, 3)) = R1 * 1.0; P3.col(3) = -t * 1.0;
        P4(Range::all(), Range(0, 3)) = R2 * 1.0; P4.col(3) = -t * 1.0;

        // Do the cheirality check.
        // Notice here a threshold dist is used to filter
        // out far away points (i.e. infinite points) since
        // there depth may vary between postive and negtive.
        double dist = 50.0;
        Mat Q;
        triangulatePoints(P0, P1, points1, points2, Q);
        Mat mask1 = Q.row(2).mul(Q.row(3)) > 0;
        Q.row(0) /= Q.row(3);
        Q.row(1) /= Q.row(3);
        Q.row(2) /= Q.row(3);
        Q.row(3) /= Q.row(3);
        mask1 = (Q.row(2) < dist) & mask1;
        Q = P1 * Q;
        mask1 = (Q.row(2) > 0) & mask1;
        mask1 = (Q.row(2) < dist) & mask1;

        triangulatePoints(P0, P2, points1, points2, Q);
        Mat mask2 = Q.row(2).mul(Q.row(3)) > 0;
        Q.row(0) /= Q.row(3);
        Q.row(1) /= Q.row(3);
        Q.row(2) /= Q.row(3);
        Q.row(3) /= Q.row(3);
        mask2 = (Q.row(2) < dist) & mask2;
        Q = P2 * Q;
        mask2 = (Q.row(2) > 0) & mask2;
        mask2 = (Q.row(2) < dist) & mask2;

        triangulatePoints(P0, P3, points1, points2, Q);
        Mat mask3 = Q.row(2).mul(Q.row(3)) > 0;
        Q.row(0) /= Q.row(3);
        Q.row(1) /= Q.row(3);
        Q.row(2) /= Q.row(3);
        Q.row(3) /= Q.row(3);
        mask3 = (Q.row(2) < dist) & mask3;
        Q = P3 * Q;
        mask3 = (Q.row(2) > 0) & mask3;
        mask3 = (Q.row(2) < dist) & mask3;

        triangulatePoints(P0, P4, points1, points2, Q);
        Mat mask4 = Q.row(2).mul(Q.row(3)) > 0;
        Q.row(0) /= Q.row(3);
        Q.row(1) /= Q.row(3);
        Q.row(2) /= Q.row(3);
        Q.row(3) /= Q.row(3);
        mask4 = (Q.row(2) < dist) & mask4;
        Q = P4 * Q;
        mask4 = (Q.row(2) > 0) & mask4;
        mask4 = (Q.row(2) < dist) & mask4;

        mask1 = mask1.t();
        mask2 = mask2.t();
        mask3 = mask3.t();
        mask4 = mask4.t();

        // If _mask is given, then use it to filter outliers.
        if (!_mask.empty())
        {
            Mat mask = _mask.getMat();
            CV_Assert(mask.size() == mask1.size());
            bitwise_and(mask, mask1, mask1);
            bitwise_and(mask, mask2, mask2);
            bitwise_and(mask, mask3, mask3);
            bitwise_and(mask, mask4, mask4);
        }
        if (_mask.empty() && _mask.needed())
        {
            _mask.create(mask1.size(), CV_8U);
        }

        CV_Assert(_R.needed() && _t.needed());
        _R.create(3, 3, R1.type());
        _t.create(3, 1, t.type());

        int good1 = countNonZero(mask1);
        int good2 = countNonZero(mask2);
        int good3 = countNonZero(mask3);
        int good4 = countNonZero(mask4);

        if (good1 >= good2 && good1 >= good3 && good1 >= good4)
        {
            R1.copyTo(_R);
            t.copyTo(_t);
            if (_mask.needed()) mask1.copyTo(_mask);
            return good1;
        }
        else if (good2 >= good1 && good2 >= good3 && good2 >= good4)
        {
            R2.copyTo(_R);
            t.copyTo(_t);
            if (_mask.needed()) mask2.copyTo(_mask);
            return good2;
        }
        else if (good3 >= good1 && good3 >= good2 && good3 >= good4)
        {
            t = -t;
            R1.copyTo(_R);
            t.copyTo(_t);
            if (_mask.needed()) mask3.copyTo(_mask);
            return good3;
        }
        else
        {
            t = -t;
            R2.copyTo(_R);
            t.copyTo(_t);
            if (_mask.needed()) mask4.copyTo(_mask);
            return good4;
        }
    }

    //简化版的recoverPose，只提供焦距和主点坐标。简化了相机模型。
    int recoverPose( InputArray E, InputArray _points1, InputArray _points2, OutputArray _R,
                         OutputArray _t, double focal, Point2d pp, InputOutputArray _mask)
    {
        Mat cameraMatrix = (Mat_<double>(3,3) << focal, 0, pp.x, 0, focal, pp.y, 0, 0, 1);
        return cv::recoverPose(E, _points1, _points2, cameraMatrix, _R, _t, _mask);
    }
}

//一共定义了两个函数，分别是decomposeEssentialMat——用于分解本质矩阵，另一个是recoverPose，用于恢复相机位姿。
//recoverPose有两个版本，第一个是需要的相机内参矩阵，第二个是简化版本，简化了输入，但是依然是调用第一个版本的逻辑
```



###### solveRelativeRT

基于特征点的两帧相对位姿（旋转和平移）求解函数。

```cpp
bool MotionEstimator::solveRelativeRT(const vector<pair<Vector3d, Vector3d>> &corres, Matrix3d &Rotation, Vector3d &Translation)
{
    if (corres.size() >= 15)
    {
        vector<cv::Point2f> ll, rr;
        //将eigen3D点转化为opencv的2D点
        for (int i = 0; i < int(corres.size()); i++)
        {
            ll.push_back(cv::Point2f(corres[i].first(0), corres[i].first(1)));
            rr.push_back(cv::Point2f(corres[i].second(0), corres[i].second(1)));
        }
        //利用RANSAC一致性算法，随机从两帧图像的点队中选取组成基础函数，然后求重投影误差，多次重复找出内点数量最多对应的基础矩阵F。
        cv::Mat mask;
        cv::Mat E = cv::findFundamentalMat(ll, rr, cv::FM_RANSAC, 0.3 / 460, 0.99, mask);
        cv::Mat cameraMatrix = (cv::Mat_<double>(3, 3) << 1, 0, 0, 0, 1, 0, 0, 0, 1);//构造3*3的单位相机内参矩阵
        cv::Mat rot, trans;
        int inlier_cnt = cv::recoverPose(E, ll, rr, cameraMatrix, rot, trans, mask);//从匹配点中恢复相机的相对位姿
        //cout << "inlier_cnt " << inlier_cnt << endl;

        Eigen::Matrix3d R;
        Eigen::Vector3d T;
        //将opencv数据转化为eigen格式
        for (int i = 0; i < 3; i++)
        {   
            T(i) = trans.at<double>(i, 0);
            for (int j = 0; j < 3; j++)
                R(i, j) = rot.at<double>(i, j);
        }

        Rotation = R.transpose();
        //前一帧 -> 当前帧 转化为 当前帧 -> 前一帧
        Translation = -R.transpose() * T;
        if(inlier_cnt > 12)
            return true;
        else
            return false;
    }
    return false;
}
```




