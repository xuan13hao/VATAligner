#ifndef __ALIGNMODEL_H__
#define __ALIGNMODEL_H__

#include "CreateDB.h"
#include "DNAAlignFlow.h"
#include "ProteinAlignFlow.h"


class RunModel
{
    public:
    void static CreateDNADB()
    {
        CreateDB(DNA());
    }
    /*
    void static CreateDoubleStrandDNADB()
    {
        CreateDNA_DB();
    }
    */
    void static CreateProteinDB()
    {
        CreateDB(Protein());
    }

    void static DNAAlign()
    {
        DNAMasterThread<DNA>();
    }

    void static ProteinAlign()
    {
        ProteinMasterThread<Protein>();
    }
};


#endif // __ALIGNMODEL_H__