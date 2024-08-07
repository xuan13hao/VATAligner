#ifndef __ALIGNMODEL_H__
#define __ALIGNMODEL_H__

#include "CreateDB.h"
#include "DNAAlignFlow.h"
#include "ProteinAlignFlow.h"
#include "BlastxAlignFlow.h"

class RunModel
{
    public:
    void static CreateDNADB()
    {
        CreateDB(DNA());
    }
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
    void static BlastxAlign()
    {
        BlastxMasterThread<Protein>();
    }
};


#endif // __ALIGNMODEL_H__