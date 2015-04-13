#pragma once

#include "ig_structs/ig_structure_structs.hpp"
#include "ig_structs/shm_settings.hpp"
#include "ig_structs/cdr_settings.hpp"

template<class VDJ_Recombination_Ptr>
class IgVariableRegion {
    VDJ_Recombination_Ptr vdj_recombination_;
    SHMSettings shm_settings_;
    CDRSettings cdr_settings_;

    string sequence_;
    bool update_seq_;

    void ComputeSequence() {
        // do something
        update_seq_ = false;
    }

public:
    IgVariableRegion(VDJ_Recombination_Ptr vdj_recombination) :
            vdj_recombination_(vdj_recombination),
            sequence_(vdj_recombination->Sequence()),
            update_seq_(false) {
    }

    VDJ_Recombination_Ptr VDJ_Recombination() {
        return vdj_recombination_;
    }

    void SetSHMSettings(SHMSettings shm_settings) {
        shm_settings_ = shm_settings;
        update_seq_ = true;
    }

    const SHMSettings& GetSHMSettings() const { return shm_settings_; };

    void SetCDRSettings(CDRSettings cdr_settings) {
        cdr_settings_ = cdr_settings;
    }

    const CDRSettings& GetCDRSettings() const { return cdr_settings_; }

    string Sequence() {
        if(update_seq_)
            ComputeSequence();
        return sequence_;
    }

    typedef shared_ptr<IgVariableRegion<VDJ_Recombination_Ptr> > IgVariableRegionPtr;

    IgVariableRegionPtr Clone() {
        IgVariableRegionPtr variable_region_ptr(
                new IgVariableRegion<VDJ_Recombination_Ptr>(vdj_recombination_->Clone()));
        variable_region_ptr->SetSHMSettings(shm_settings_);
        variable_region_ptr->SetCDRSettings(cdr_settings_);
        return variable_region_ptr;
    }
};

typedef IgVariableRegion<HC_VDJ_Recombination_Ptr> HC_VariableRegion;
typedef IgVariableRegion<LC_VDJ_Recombination_Ptr> LC_VariableRegion;
typedef shared_ptr<IgVariableRegion<HC_VDJ_Recombination_Ptr> > HC_VariableRegionPtr;
typedef shared_ptr<IgVariableRegion<LC_VDJ_Recombination_Ptr> > LC_VariableRegionPtr;

// ----------------------------------------------------------------------------

// cluster is characterized by:
//      - IgVariableRegion (VDJ recombination, CDR3 labeling)
//      - multiplicity

template<class IgVariableRegionPtr>
class IgCluster {
    IgVariableRegionPtr variable_region_;
    size_t multiplicity_;

public:
    IgCluster(IgVariableRegionPtr variable_region, size_t multiplicity) :
            variable_region_(variable_region),
            multiplicity_(multiplicity) { }

    IgVariableRegionPtr IgVariableRegion() const { return variable_region_; }

    size_t Multiplicity() const { return multiplicity_; }
};

typedef IgCluster<HC_VariableRegionPtr> HC_Cluster;
typedef IgCluster<LC_VariableRegionPtr> LC_Cluster;

typedef vector<HC_Cluster>::const_iterator HC_ClusterIterator;
typedef vector<LC_Cluster>::const_iterator LC_ClusterIterator;

// ----------------------------------------------------------------------------

template<class IgCluster, class IgClusterIterator>
class Repertoire {
    vector<IgCluster> ig_clusters_;

public:
    Repertoire() { }

    void Add(IgCluster ig_cluster, size_t cluster_multiplicity = 1) {
        for(size_t i = 0; i < cluster_multiplicity; i++)
            ig_clusters_.push_back(ig_cluster);
    }

    size_t Size() const {
        return ig_clusters_.size();
    }

    IgClusterIterator begin() const {
        return ig_clusters_.begin();
    }

    IgClusterIterator end() const {
        return ig_clusters_.end();
    }
};

typedef Repertoire<HC_Cluster, HC_ClusterIterator> HC_Repertoire;
typedef Repertoire<LC_Cluster, LC_ClusterIterator> LC_Repertoire;