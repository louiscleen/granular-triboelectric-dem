#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include "Disk.h"
#include "Cell.h"

using namespace std;

Disk::Disk()
{
    m_index = -1;
    m_linkedDisk = nullptr;
    m_linkedCell = nullptr;
    m_radius = 0.;
    m_mass = 0.;
    m_inertia = 0.5*m_mass*m_radius*m_radius;
    
    m_r = Eigen::Vector3d::Zero();
    m_v = Eigen::Vector3d::Zero();
    m_a = Eigen::Vector3d::Zero();
    m_F = Eigen::Vector3d::Zero();
    
    m_theta = 0.;
    m_w = Eigen::Vector3d::Zero();
    m_alpha = Eigen::Vector3d::Zero();
    m_M = Eigen::Vector3d::Zero();

    m_patch = 0; // Nombre de patchs
    m_sat = 0; // Nombre maximum de charges par patch
    m_patch_state = vector<bool>(m_patch,false);
    m_patch_charges = vector<int>(m_patch,0);
    //std::vector<double> m_patch_angle(i_p);
}

Disk::Disk(int i_index, int i_patch, int i_sat, vector<bool> i_patch_state, double i_radius, double i_mass, double i_x, double i_y, double i_vx, double i_vy)
{
    m_index = i_index;
    m_linkedDisk = nullptr;
    m_linkedCell = nullptr;
    m_radius = i_radius;
    m_mass = i_mass;
    m_inertia = 0.5*m_mass*m_radius*m_radius;
    
    m_r = Eigen::Vector3d(i_x,i_y,0.);
    m_v = Eigen::Vector3d(i_vx,i_vy,0.);
    m_a = Eigen::Vector3d::Zero();
    m_F = Eigen::Vector3d::Zero();
    
    m_theta = 0.;
    m_w = Eigen::Vector3d::Zero();
    m_alpha = Eigen::Vector3d::Zero();
    m_M = Eigen::Vector3d::Zero();

    m_patch = i_patch; // Nombre de patchs
    m_sat = i_sat; // Nombre maximum de charges par patch
    m_patch_state = i_patch_state;
    m_patch_charges = vector<int>(i_patch,0);
    //std::vector<double> m_patch_angle(i_p);
}

Disk::~Disk()
{
    
}

void Disk::reset_force()
{
    m_F = Eigen::Vector3d::Zero();
    m_M = Eigen::Vector3d::Zero();
}

void Disk::add_force(const Eigen::Vector3d& i_F)
{
    m_F += i_F;
}

void Disk::add_momentum(const Eigen::Vector3d& i_M)
{
    m_M += i_M;
}

void Disk::add_gravity_force(const Eigen::Vector3d& i_g)
{
    m_F += m_mass*i_g;
}

void Disk::update_position(double dt)
{
    m_r += m_v*dt;
    m_theta += m_w.z()*dt;
}

void Disk::update_velocity(double dt)
{
    m_a = m_F/m_mass;
    m_v += m_a*dt;
    
    m_alpha = m_M/m_inertia;
    m_w += m_alpha*dt;
}

void Disk::set_linked_disk(Disk* i_linkedDisk)
{
    m_linkedDisk = i_linkedDisk;
}

void Disk::set_linked_cell(Cell* i_linkedCell)
{
    m_linkedCell = i_linkedCell;
}


void Disk::print(int i_num)
{
    namespace fs = std::filesystem;

    // Vérifie si le dossier 'data' existe, sinon le crée
    if (!fs::exists("data")) {
        fs::create_directory("data");
    }

    ofstream myfile;
    string fileName = "data/grain" + to_string(i_num) + ".txt";
    myfile.open(fileName, ios::app);
    myfile.precision(10);
    myfile << m_index << "\t" 
           << m_r.x() << "\t" << m_r.y() << "\t" 
           << m_v.x() << "\t" << m_v.y() << "\t" 
           << m_theta << "\t" << m_radius << endl;
    myfile.close();
}



Disk* Disk::linked_disk()
{
    return m_linkedDisk;
}

Cell* Disk::linked_cell()
{
    return m_linkedCell;
}

int Disk::index()
{
    return m_index;
}

double Disk::radius()
{
    return m_radius;
}

double Disk::mass()
{
    return m_mass;
}

Eigen::Vector3d Disk::r()
{
    return m_r;
}

Eigen::Vector3d Disk::v()
{
    return m_v;
}

double Disk::theta()
{
    return m_theta;
}

Eigen::Vector3d Disk::w()
{
    return m_w;
}

bool Disk::is_touching(double i_x,double i_y, double i_radius)
{
    return (m_r - Eigen::Vector3d(i_x,i_y,0.)).norm() < (m_radius+i_radius);
}

void Disk::set_position(double x, double y)
{
    m_r.x() = x; // Correctly set the x component
    m_r.y() = y; // Correctly set the y component
}

void Disk::setCaged(bool state) 
{ 
    caged = state; 
}

bool Disk::isCaged() 
{ 
    return caged; 
}

int Disk::patch_count() 
{ 
    return m_patch; 
}

std::vector<bool> Disk::getPatchStates() 
{ 
    return m_patch_state; 
}

void Disk::set_charge(int patch_index, int charge) 
{ 
    m_patch_charges[patch_index] = charge; 
}

void Disk::add_charge(int patch_index, int charge) 
{ 
    m_patch_charges[patch_index] += charge; 
}

void Disk::reset_patch_charges() 
{ 
    std::fill(m_patch_charges.begin(), m_patch_charges.end(), 0); 
}

void Disk::reset_patch_charge(int patch_index) 
{ 
    m_patch_charges[patch_index] = 0; 
}

std::vector<int> Disk::getPatchCharges() 
{ 
    return m_patch_charges; 
}

int Disk::getPatchCharge(int patch_index) 
{ 
    return m_patch_charges[patch_index]; 
}

int Disk::get_sat() 
{ 
    return m_sat; 
}


void Disk::addContact(int myPatchIndex, int neighborIndex, int neighborPatchIndex, bool state) 
{
    // Vérifie si le contact existe déjà
    for (auto& contact : contacts) 
    {
        if (contact.patchIndex == myPatchIndex && contact.neighborIndex == neighborIndex && contact.neighborPatchIndex == neighborPatchIndex) {
            contact.state = state; // Met à jour l'état du contact
            return; // Le contact existe déjà, ne rien faire
        }
    }
    
    // Ajoute un nouveau contact
    ContactInfo newContact;
    newContact.patchIndex = myPatchIndex;
    newContact.neighborIndex = neighborIndex;
    newContact.neighborPatchIndex = neighborPatchIndex;
    newContact.state = state;
    newContact.hasTransferredCharge = false;
    contacts.push_back(newContact);
}


std::vector<Disk::ContactInfo>& Disk::getContacts() 
{
    return contacts;
}

bool Disk::getContactState(int myPatchIndex, int neighborIndex, int neighborPatchIndex) 
{
    for (const auto& contact : contacts) 
    {
        if (contact.patchIndex == myPatchIndex && contact.neighborIndex == neighborIndex && contact.neighborPatchIndex == neighborPatchIndex) {
            return contact.state; // Retourne l'état du contact
        }
    }
    return false; // Retourne false si le contact n'existe pas
}

void Disk::setContactState(int myPatchIndex, int neighborIndex, int neighborPatchIndex, bool state) 
{
    for (auto& contact : contacts) 
    {
        if (contact.patchIndex == myPatchIndex && contact.neighborIndex == neighborIndex && contact.neighborPatchIndex == neighborPatchIndex) {
            contact.state = state; // Met à jour l'état du contact
            return;
        }
    }
}

bool Disk::getContactChargeTransfer(int myPatchIndex, int neighborIndex, int neighborPatchIndex) 
{
    for (const auto& contact : contacts) 
    {
        if (contact.patchIndex == myPatchIndex && contact.neighborIndex == neighborIndex && contact.neighborPatchIndex == neighborPatchIndex) {
            return contact.hasTransferredCharge; // Retourne l'état du transfert de charge
        }
    }
    return false; // Retourne false si le contact n'existe pas
}

void Disk::setContactChargeTransfer(int myPatchIndex, int neighborIndex, int neighborPatchIndex, bool hasTransferredCharge) 
{
    for (auto& contact : contacts) 
    {
        if (contact.patchIndex == myPatchIndex && contact.neighborIndex == neighborIndex && contact.neighborPatchIndex == neighborPatchIndex) {
            contact.hasTransferredCharge = hasTransferredCharge; // Met à jour l'état du transfert de charge
            return;
        }
    }
}

double Disk::inertia()
{
    return m_inertia;
}