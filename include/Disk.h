#ifndef DISK_H
#define DISK_H

#include <Eigen/Dense>
#include <vector>

class Cell;

class Disk
{
private:
    int m_index;
    Cell* m_linkedCell;
    Disk* m_linkedDisk;
    
    double m_radius, m_mass, m_inertia;
    
    //translation
    Eigen::Vector3d m_r;
    Eigen::Vector3d m_v;
    Eigen::Vector3d m_a;
    Eigen::Vector3d m_F;

    //rotation
    double m_theta;
    Eigen::Vector3d m_w;
    Eigen::Vector3d m_alpha;
    Eigen::Vector3d m_M;

    bool caged;
    int m_patch; // Nombre de patchs
    int m_sat; // Nombre maximum de charges par patch
    

    
    std::vector<bool> m_patch_state; // Tableau des états des patchs (1 = accepteur, 0 = donneur)
    std::vector<int> m_patch_charges; // Tableau des charges des patchs



public:
    Disk();
    Disk(int, int, int, std::vector<bool>, double, double, double, double, double, double);
    ~Disk();
    
    void update_velocity(double);
    void update_position(double);
    void reset_force();
    void add_force(const Eigen::Vector3d&);
    void add_momentum(const Eigen::Vector3d&);
    void add_gravity_force(const Eigen::Vector3d&);
    void set_linked_disk(Disk*);
    void set_linked_cell(Cell*);
    void print(int);
    void set_position(double x, double y);
    
    Cell* linked_cell();
    Disk* linked_disk();
    
    int index();
    double radius();
    double mass();
    double theta();
    Eigen::Vector3d r();
    Eigen::Vector3d v();
    Eigen::Vector3d w();

    struct ContactInfo {
        int patchIndex;               // Index du patch en contact
        int neighborIndex;            // Index de la particule en contact
        int neighborPatchIndex;       // Index du patch de la particule en contact
        bool state;                   // Si le transfert est récent
        bool hasTransferredCharge;    // Si un transfert de charge a eu lieu

    };

    std::vector<ContactInfo> contacts;
    
    bool is_touching(double,double,double);

    void setCaged(bool state);
    bool isCaged();

    int patch_count();
    std::vector<bool> getPatchStates();
    void add_charge(int, int);
    std::vector<int> getPatchCharges();
    int getPatchCharge(int);
    int get_sat();
    

    void addContact(int, int, int, bool);
    std::vector<ContactInfo>& getContacts();
    bool getContactState(int, int, int);
    void setContactState(int, int, int, bool);
    bool getContactChargeTransfer(int, int, int);
    void setContactChargeTransfer(int, int, int, bool);
    double inertia();

};

#endif // DISK_H