#include "Plan.h"
#include <math.h>

Plan::Plan(double i_x, double i_y, double i_nx, double i_ny) {
    m_r = Eigen::Vector3d(i_x, i_y, 0.);
    m_n = Eigen::Vector3d(i_nx, i_ny, 0.);
}

Plan::~Plan() {}

Eigen::Vector3d Plan::r() { return m_r; }

Eigen::Vector3d Plan::n() { return m_n; }

void Plan::set_position(double x, double y) {
    m_r = Eigen::Vector3d(x, y, 0.); // ou selon la déclaration exacte de r_
}