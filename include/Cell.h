#ifndef CELL_H
#define CELL_H

#include <vector>

class Disk;

class Cell {

  private:
    int m_index;
    Disk *m_head_of_list;
    std::vector<Cell *> m_neighbors;
    std::vector<Cell *> m_neighborsOfNeighbors;

  public:
    Cell(int);
    ~Cell();

    void set_head_of_list(Disk *);
    void add_neighbor(Cell *);
    void add_neighbor_of_neighbor(Cell *);

    Disk *head_of_list();
    std::vector<Cell *> neighbors();
    std::vector<Cell *> neighbors_of_neighbors();
};

#endif // CELL_H