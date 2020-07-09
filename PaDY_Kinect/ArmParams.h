#pragma once
#define _USE_MATH_DEFINES
#include<math.h>

#ifndef DEG2RAD
#define DEG2RAD    (M_PI/180.0)
#endif

#ifndef Deg2Rad
#define Deg2Rad    (M_PI/180.0)
#endif

#ifndef RAD2DEG
#define RAD2DEG    (180.0/M_PI)
#endif

#ifndef Rad2Deg
#define Rad2Deg    (180.0/M_PI)
#endif


#define J_NUM 2  //Total number of joints
#define M1 (5.0) //Mass of link1
#define M2 (5.5) //Mass of link2
#define L1 (1.068) //Length of link1
#define L2 (0.940) // Length of link2
#define LG1 (L1 / 2.0) //Center of link1
#define LG2 (L2 / 2.0) //Center of link2
#define INERTIA1 (M1 * L1 * L1 / 3.0) //Inertia of link1
#define INERTIA2 (M2 * L2 * L2 / 3.0) //Inertia of link2
#define ANVEL1_MAX (M_PI) //Max angular velocity of joint1
#define ANVEL2_MAX (M_PI) //Max angular velocity of joint2
#define ANACC1_MAX (M_PI) //Max angular accelaration of joint1
#define ANACC2_MAX (M_PI) //Max angular accelaration of joint2