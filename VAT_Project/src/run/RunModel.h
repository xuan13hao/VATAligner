#ifndef __ALIGNMODEL_H__
#define __ALIGNMODEL_H__

#include "CreateDB.h"
#include "DNAAlignFlow.h"
#include "ProteinAlignFlow.h"


class RunModel
{
    public:
    void CreateDNADB()
    {
        CreateDB(DNA());
    }
    void CreateProteinDB()
    {
        CreateDB(Protein());
    }

    void DNAAlign()
    {
        DNAMasterThread<DNA>();
    }

    void ProteinAlign()
    {
        ProteinMasterThread<Protein>();
    }

};


#endif // __ALIGNMODEL_H__