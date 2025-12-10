##### rosNodeTest.cpp

在学习这部分代码的时候，需要先了解一个概念**回调函数**。这是一个广泛存在于编程中的概念，常用于事件驱动编程或异步处理。在ROS中，回调函数的使用常用于消息订阅（消息到达调用回调函数）和定时器（到时间调用回调函数）。

整个源文件通过定义多个回调函数，并在主函数中创建节点和订阅话题，形成一个 VIO 数据处理模块。该节点通过回调函数接收传感器信息，再通过 estimator 接口执行实际的状态估计和计算。

```cpp
#include <stdio.h>
#include <queue>
#include <map>
#include <thread>
#include <mutex>
#include <ros/ros.h>
#include <cv_bridge/cv_bridge.h>
#include <opencv2/opencv.hpp>
#include "estimator/estimator.h"
#include "estimator/parameters.h"
#include "utility/visualization.h"

Estimator estimator;

//定义了节点中用来缓存传感器消息的队列，意思是将收到的消息先存储起来，再在同步线程里面统一处理。
//为什么要存起来？
//因为IMU、相机的数据采集频率是不一样的，要保证在同一时间处理对应的时间戳的数据，则要先将数据存储起来。
queue<sensor_msgs::ImuConstPtr> imu_buf;//缓存IMU消息
queue<sensor_msgs::PointCloudConstPtr> feature_buf;//缓存视觉特征点消息
queue<sensor_msgs::ImageConstPtr> img0_buf;//缓存左相机图像
queue<sensor_msgs::ImageConstPtr> img1_buf;//缓存右相机图像（双目使用）
std::mutex m_buf;//互斥锁，保证多线程同时访问队列安全
```

###### img0_callback

ROS回调函数，当节点收到左相机图像时，把图像消息存入队列。

```cpp
void img0_callback(const sensor_msgs::ImageConstPtr &img_msg)
{
    m_buf.lock();//数据上锁
    img0_buf.push(img_msg);//将图像存入队列
    m_buf.unlock();//数据解锁
}
```

###### img1_callback

同上，处理右图像。

```cpp
void img1_callback(const sensor_msgs::ImageConstPtr &img_msg)
{
    m_buf.lock();
    img1_buf.push(img_msg);
    m_buf.unlock();
}
```

###### getImageFromMsg

该函数将ROS图像消息转化为OpenCV的 cv::Mat格式。

sensor_msgs::ImageConstPtr 是ROS里图像消息的智能指针类型。

```cpp
cv::Mat getImageFromMsg(const sensor_msgs::ImageConstPtr &img_msg)
{
    cv_bridge::CvImageConstPtr ptr;
    //对于8UC1类型，cv_brigde内部没有对应的解析规则，所以需要手动的复制信息。对于其他类型，则能够自动复制信息。
    if (img_msg->encoding == "8UC1")
    {
        sensor_msgs::Image img;
        img.header = img_msg->header;
        img.height = img_msg->height;
        img.width = img_msg->width;
        img.is_bigendian = img_msg->is_bigendian;
        img.step = img_msg->step;
        img.data = img_msg->data;
        img.encoding = "mono8";
        ptr = cv_bridge::toCvCopy(img, sensor_msgs::image_encodings::MONO8);
    }
    else
        ptr = cv_bridge::toCvCopy(img_msg, sensor_msgs::image_encodings::MONO8);

    cv::Mat img = ptr->image.clone();
    return img;
}
```

###### sync_process

同步线程函数，从缓存队列中取出时间对齐的图像并送给Estimator处理。

```cpp
// extract images with same timestamp from two topics
void sync_process()
{
    while(1)
    {
        if(STEREO)//双目同步处理
        {
            cv::Mat image0, image1;
            std_msgs::Header header;
            //std_msgs是ROS提供的一个标准消息库，定义了一些常用的消息类型。
            //Header消息头，包含时间戳，序号和坐标系
            double time = 0;
            m_buf.lock();
            if (!img0_buf.empty() && !img1_buf.empty())
            {
                double time0 = img0_buf.front()->header.stamp.toSec();
                //获取左图像的时间戳
                double time1 = img1_buf.front()->header.stamp.toSec();
                //获取右图像时间戳
                
                // 0.003s sync tolerance
                //时间容差是0.003S
                if(time0 < time1 - 0.003)
                {
                    img0_buf.pop();
                    printf("throw img0\n");
                }
                else if(time0 > time1 + 0.003)
                {
                    img1_buf.pop();
                    printf("throw img1\n");
                }
                else
                {
                    //左右图像配对成功之后，获取时间戳，并将左右图转化为矩阵形式存储，然后从队列中删除。
                    time = img0_buf.front()->header.stamp.toSec();
                    header = img0_buf.front()->header;
                    image0 = getImageFromMsg(img0_buf.front());
                    img0_buf.pop();
                    image1 = getImageFromMsg(img1_buf.front());
                    img1_buf.pop();
                    //printf("find img0 and img1\n");
                }
            }
            m_buf.unlock();
            //确保左右图不为空后送入Estimator
            if(!image0.empty())
                estimator.inputImage(time, image0, image1);
        }
        else
            //处理单目相机
        {
            cv::Mat image;
            std_msgs::Header header;
            double time = 0;
            m_buf.lock();
            if(!img0_buf.empty())
            {
                time = img0_buf.front()->header.stamp.toSec();
                header = img0_buf.front()->header;
                image = getImageFromMsg(img0_buf.front());
                img0_buf.pop();
            }
            m_buf.unlock();
            if(!image.empty())
                estimator.inputImage(time, image);
        }

        //让同步线程休息2ms，避免无线循环占用过多CPU
        std::chrono::milliseconds dura(2);//定义睡眠时间
        std::this_thread::sleep_for(dura);//线程休眠
    }
}
```

###### imu_callback

IMU消息回调函数，负责接受ROS发布的IMU数据并送入Estimator。

sensor_msgs是ROS的标准传感器消息库，提供了一系列常用的传感器数据类型。其包含的消息类型有Imu、Image、CameraInfo、laserScan、PointCloud(2)、Range、JointState等。其中IMU数据包含了加速度、角速度和姿态四元数。

```cpp
void imu_callback(const sensor_msgs::ImuConstPtr &imu_msg)
{
    double t = imu_msg->header.stamp.toSec();
    double dx = imu_msg->linear_acceleration.x;
    double dy = imu_msg->linear_acceleration.y;
    double dz = imu_msg->linear_acceleration.z;
    double rx = imu_msg->angular_velocity.x;
    double ry = imu_msg->angular_velocity.y;
    double rz = imu_msg->angular_velocity.z;
    Vector3d acc(dx, dy, dz);
    Vector3d gyr(rx, ry, rz);
    estimator.inputIMU(t, acc, gyr);
    return;
}
```

###### feature_callback

特征点消息回调函数，用于接受视觉跟踪器发布的特征点信息，并且处理成Estimator能够处理的格式。

```cpp
void feature_callback(const sensor_msgs::PointCloudConstPtr &feature_msg)
{
    map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>> featureFrame;
    for (unsigned int i = 0; i < feature_msg->points.size(); i++)
    {
        int feature_id = feature_msg->channels[0].values[i];
        int camera_id = feature_msg->channels[1].values[i];
        double x = feature_msg->points[i].x;
        double y = feature_msg->points[i].y;
        double z = feature_msg->points[i].z;
        double p_u = feature_msg->channels[2].values[i];
        double p_v = feature_msg->channels[3].values[i];
        double velocity_x = feature_msg->channels[4].values[i];
        double velocity_y = feature_msg->channels[5].values[i];
        //如果通道数大于5，说明存储了真值，常用于进行精度评估。
        if(feature_msg->channels.size() > 5)
        {
            double gx = feature_msg->channels[6].values[i];
            double gy = feature_msg->channels[7].values[i];
            double gz = feature_msg->channels[8].values[i];
            pts_gt[feature_id] = Eigen::Vector3d(gx, gy, gz);
            //printf("receive pts gt %d %f %f %f\n", feature_id, gx, gy, gz);
            
        }
        ROS_ASSERT(z == 1);//确保归一化的平面坐标z==1
        Eigen::Matrix<double, 7, 1> xyz_uv_velocity;
        xyz_uv_velocity << x, y, z, p_u, p_v, velocity_x, velocity_y;
        featureFrame[feature_id].emplace_back(camera_id,  xyz_uv_velocity);
    }
    double t = feature_msg->header.stamp.toSec();
    estimator.inputFeature(t, featureFrame);
    return;
}

```

###### restart_callback

重启回调函数，在接收到重启命令时，重新初始化估计器。

```cpp
void restart_callback(const std_msgs::BoolConstPtr &restart_msg)
{
    if (restart_msg->data == true)
    {
        ROS_WARN("restart the estimator!");
        estimator.clearState();
        estimator.setParameter();
    }
    return;
}
```

###### imu_switch_callback

IMU开关回调函数，用于在运行过程中根据收到的消息开启或者关闭IMU传感器的数据融合。

```cpp
void imu_switch_callback(const std_msgs::BoolConstPtr &switch_msg)
{
    if (switch_msg->data == true)
    {
        //ROS_WARN("use IMU!");
        estimator.changeSensorType(1, STEREO);
    }
    else
    {
        //ROS_WARN("disable IMU!");
        estimator.changeSensorType(0, STEREO);
    }
    return;
}
```

###### cam_switch_callback

相机单双目切换回调函数

```cpp
void cam_switch_callback(const std_msgs::BoolConstPtr &switch_msg)
{
    if (switch_msg->data == true)
    {
        //ROS_WARN("use stereo!");
        estimator.changeSensorType(USE_IMU, 1);
    }
    else
    {
        //ROS_WARN("use mono camera (left)!");
        estimator.changeSensorType(USE_IMU, 0);
    }
    return;
}
```

###### 主函数

负责初始化ROS节点、加载函数、订阅话题、启动同步线程，并进入主循环。

什么是节点？

节点（Node）是ROS系统的最小运行单位。负责处理某一功能或者任务。可以用于数据获取、处理、发布、服务调用与提供等。

节点本质上是一个独立运行的程序，复杂执行特定功能。一个节点就像一个“功能模块”，强调独立功能与模块化，在C++中没有这个概念，本身是分布式系统中的一个概念，在ROS中用于实现模块化和通信。

类似的在ROS中有一些概念，下面通过C++里面的一些概念对ROS里面的概念进行类别，帮助理解。

|            ROS概念             |           C++类比            |
| :----------------------------: | :--------------------------: |
|          节点（Node）          | 独立的类或可执行程序（.exe） |
|     节点句柄（NodeHandle)      |     指向对象的指针或引用     |
|        消息（Message）         |        结构体/类对象         |
|      发布者（Publisher）       |   写入数据的函数或对象方法   |
|      订阅者（Subscriber）      |   读取数据的函数或回调函数   |
|        服务（Service）         |           函数调用           |
|         话题（Topic）          |           全局变量           |
| 参数服务器（Parameter Server） |      全局变量或配置文件      |

服务和话题揭示了ROS两种主要的通信方式，话题通信是异步流式通信（单向的消息传递，不需要回应），服务通信是同步请求-响应通信（双向的，需要反馈）。

```cpp
int main(int argc, char **argv)
{
    ros::init(argc, argv, "vins_estimator");//初始化ROS节点。
    ros::NodeHandle n("~");//创建节点句柄，传入“~”表示私有命名空间。
    //节点句柄（NodeHandle）是节点访问ROS系统的入口
    ros::console::set_logger_level(ROSCONSOLE_DEFAULT_NAME, ros::console::levels::Info);//设置当前节点的日志输出等级

    //提示用户必须要有一个配置文件路径才能够运行
    if(argc != 2)
    {
        printf("please intput: rosrun vins vins_node [config file] \n"
               "for example: rosrun vins vins_node "
               "~/catkin_ws/src/VINS-Fusion/config/euroc/euroc_stereo_imu_config.yaml \n");
        return 1;
    }

    //从命令行参数中读取配置文件路径
    string config_file = argv[1];
    printf("config_file: %s\n", argv[1]);

    //从路径中读取参数，并进行存储
    readParameters(config_file);
    estimator.setParameter();

    //禁止Eigen多线程运算
#ifdef EIGEN_DONT_PARALLELIZE
    ROS_DEBUG("EIGEN_DONT_PARALLELIZE");
#endif

    ROS_WARN("waiting for image and imu...");

    //复杂创建需要发布的ROS输出话题。
    //简单的理解为将计算得到的参数都集中起来，方便传递。
    registerPub(n);

    //定义一个ROS订阅者对象，用于接收IMU传感器数据
    ros::Subscriber sub_imu;
    if(USE_IMU)
    {
        sub_imu = n.subscribe(IMU_TOPIC, 2000, imu_callback, ros::TransportHints().tcpNoDelay());
    }
    ros::Subscriber sub_feature = n.subscribe("/feature_tracker/feature", 2000, feature_callback);
    ros::Subscriber sub_img0 = n.subscribe(IMAGE0_TOPIC, 100, img0_callback);
    ros::Subscriber sub_img1;
    if(STEREO)
    {
        sub_img1 = n.subscribe(IMAGE1_TOPIC, 100, img1_callback);
    }
    ros::Subscriber sub_restart = n.subscribe("/vins_restart", 100, restart_callback);
    ros::Subscriber sub_imu_switch = n.subscribe("/vins_imu_switch", 100, imu_switch_callback);
    ros::Subscriber sub_cam_switch = n.subscribe("/vins_cam_switch", 100, cam_switch_callback);

    //创建新线程进行图像同步的处理。实现了多线程并行，一个处理消息回调，一个处理图像同步。
    std::thread sync_thread{sync_process};
    ros::spin();

    return 0;
}
```

