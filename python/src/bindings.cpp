#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "llka_main.h"
#include "llka_structure.h"
#include "llka_classification.h"
#include "llka_minicif.h"
#include "llka_resource_loaders.h"
#include "llka_measurements.h"
#include "llka_nucleotide.h"
#include "llka_ntc.h"
#include "llka_superposition.h"
#include "llka_util.h"
#include "llka_connectivity_similarity.h"

namespace py = pybind11;

// Exception class for LLKA errors
class LLKAError : public std::runtime_error {
public:
    explicit LLKAError(const std::string& msg) : std::runtime_error(msg) {}
};

// Helper to check return codes and throw exceptions
void check_retcode(LLKA_RetCode code, const char* operation = nullptr) {
    if (code == LLKA_OK) return;

    std::string msg;
    if (operation) {
        msg = std::string(operation) + " failed: ";
    }

    switch (code) {
        case LLKA_E_INVALID_ARGUMENT:
            throw LLKAError(msg + "Invalid argument");
        case LLKA_E_MISSING_ATOMS:
            throw LLKAError(msg + "Missing atoms");
        case LLKA_E_MISMATCHING_SIZES:
            throw LLKAError(msg + "Mismatching sizes");
        case LLKA_E_NOT_IMPLEMENTED:
            throw LLKAError(msg + "Not implemented");
        case LLKA_E_MULTIPLE_ALT_IDS:
            throw LLKAError(msg + "Multiple alternate IDs");
        case LLKA_E_MISMATCHING_DATA:
            throw LLKAError(msg + "Mismatching data");
        case LLKA_E_BAD_CLASSIFICATION_CLUSTERS:
            throw LLKAError(msg + "Bad classification clusters");
        case LLKA_E_BAD_GOLDEN_STEPS:
            throw LLKAError(msg + "Bad golden steps");
        case LLKA_E_BAD_CONFALS:
            throw LLKAError(msg + "Bad confals");
        case LLKA_E_BAD_AVERAGE_NU_ANGLES:
            throw LLKAError(msg + "Bad average nu angles");
        case LLKA_E_BAD_CLASSIFICATION_LIMITS:
            throw LLKAError(msg + "Bad classification limits");
        case LLKA_E_NO_FILE:
            throw LLKAError(msg + "File not found");
        case LLKA_E_CANNOT_READ_FILE:
            throw LLKAError(msg + "Cannot read file");
        case LLKA_E_BAD_DATA:
            throw LLKAError(msg + "Bad data");
        case LLKA_E_NO_DATA:
            throw LLKAError(msg + "No data");
        case LLKA_E_NOTHING_TO_CLASSIFY:
            throw LLKAError(msg + "Nothing to classify");
        default:
            throw LLKAError(msg + "Unknown error (code: " + std::to_string(code) + ")");
    }
}

// RAII wrapper for LLKA_Structure
class StructureWrapper {
private:
    LLKA_Structure structure_;
    bool owns_data_;

public:
    StructureWrapper() : structure_{nullptr, 0}, owns_data_(false) {}

    explicit StructureWrapper(LLKA_Structure s) : structure_(s), owns_data_(true) {}

    // Non-owning wrapper (for borrowing from ImportedStructure)
    StructureWrapper(LLKA_Structure s, bool owns) : structure_(s), owns_data_(owns) {}

    ~StructureWrapper() {
        if (owns_data_ && structure_.atoms) {
            LLKA_destroyStructure(&structure_);
        }
    }

    // Disable copy, enable move
    StructureWrapper(const StructureWrapper&) = delete;
    StructureWrapper& operator=(const StructureWrapper&) = delete;

    StructureWrapper(StructureWrapper&& other) noexcept
        : structure_(other.structure_), owns_data_(other.owns_data_) {
        other.structure_ = {nullptr, 0};
        other.owns_data_ = false;
    }

    StructureWrapper& operator=(StructureWrapper&& other) noexcept {
        if (this != &other) {
            if (owns_data_ && structure_.atoms) {
                LLKA_destroyStructure(&structure_);
            }
            structure_ = other.structure_;
            owns_data_ = other.owns_data_;
            other.structure_ = {nullptr, 0};
            other.owns_data_ = false;
        }
        return *this;
    }

    LLKA_Structure* get() { return &structure_; }
    const LLKA_Structure* get() const { return &structure_; }

    size_t n_atoms() const { return structure_.nAtoms; }

    // Get atom by index
    const LLKA_Atom* get_atom(size_t idx) const {
        if (idx >= structure_.nAtoms) {
            throw std::out_of_range("Atom index out of range");
        }
        return &structure_.atoms[idx];
    }

    // Update coordinates from numpy array
    void update_coordinates(py::array_t<double> coords) {
        auto buf = coords.request();

        if (buf.ndim != 2) {
            throw std::runtime_error("Coordinates must be a 2D array");
        }

        if (buf.shape[0] != static_cast<ssize_t>(structure_.nAtoms)) {
            throw std::runtime_error(
                "Coordinate array must have shape (n_atoms, 3), got (" +
                std::to_string(buf.shape[0]) + ", " +
                std::to_string(buf.shape[1]) + "), expected (" +
                std::to_string(structure_.nAtoms) + ", 3)"
            );
        }

        if (buf.shape[1] != 3) {
            throw std::runtime_error("Coordinates must have 3 columns (x, y, z)");
        }

        double* ptr = static_cast<double*>(buf.ptr);

        for (size_t i = 0; i < structure_.nAtoms; ++i) {
            structure_.atoms[i].coords.x = ptr[i * 3 + 0];
            structure_.atoms[i].coords.y = ptr[i * 3 + 1];
            structure_.atoms[i].coords.z = ptr[i * 3 + 2];
        }
    }

    // Get coordinates as numpy array (copy)
    py::array_t<double> get_coordinates() const {
        py::array_t<double> result({structure_.nAtoms, size_t(3)});
        auto buf = result.request();
        double* ptr = static_cast<double*>(buf.ptr);

        for (size_t i = 0; i < structure_.nAtoms; ++i) {
            ptr[i * 3 + 0] = structure_.atoms[i].coords.x;
            ptr[i * 3 + 1] = structure_.atoms[i].coords.y;
            ptr[i * 3 + 2] = structure_.atoms[i].coords.z;
        }

        return result;
    }
};

// ImportedStructure wrapper with RAII (holds both structure and CIF data)
class ImportedStructureWrapper {
private:
    LLKA_ImportedStructure imported_;
    bool owns_data_;

public:
    ImportedStructureWrapper() : owns_data_(false) {
        memset(&imported_, 0, sizeof(imported_));
    }

    explicit ImportedStructureWrapper(LLKA_ImportedStructure imported)
        : imported_(imported), owns_data_(true) {}

    ~ImportedStructureWrapper() {
        if (owns_data_) {
            LLKA_destroyImportedStructure(&imported_);
        }
    }

    // Prevent copying
    ImportedStructureWrapper(const ImportedStructureWrapper&) = delete;
    ImportedStructureWrapper& operator=(const ImportedStructureWrapper&) = delete;

    // Allow moving
    ImportedStructureWrapper(ImportedStructureWrapper&& other) noexcept
        : imported_(other.imported_), owns_data_(other.owns_data_) {
        other.owns_data_ = false;
    }

    ImportedStructureWrapper& operator=(ImportedStructureWrapper&& other) noexcept {
        if (this != &other) {
            if (owns_data_) {
                LLKA_destroyImportedStructure(&imported_);
            }
            imported_ = other.imported_;
            owns_data_ = other.owns_data_;
            other.owns_data_ = false;
        }
        return *this;
    }

    const LLKA_ImportedStructure* get() const { return &imported_; }
    LLKA_ImportedStructure* get_mutable() { return &imported_; }

    StructureWrapper get_structure() const {
        // Return non-owning wrapper - ImportedStructure still owns the data
        return StructureWrapper(imported_.structure, false);
    }

    std::string get_entry_id() const {
        return std::string(imported_.entry.id ? imported_.entry.id : "");
    }

    LLKA_CifData* get_cif_data() const {
        return imported_.cifData;
    }
};

// RAII wrapper for ClassificationContext
class ClassificationContextWrapper {
private:
    LLKA_ClassificationContext* ctx_;

public:
    ClassificationContextWrapper() : ctx_(nullptr) {}
    explicit ClassificationContextWrapper(LLKA_ClassificationContext* ctx) : ctx_(ctx) {}

    ~ClassificationContextWrapper() {
        if (ctx_) {
            LLKA_destroyClassificationContext(ctx_);
        }
    }

    // Disable copy, enable move
    ClassificationContextWrapper(const ClassificationContextWrapper&) = delete;
    ClassificationContextWrapper& operator=(const ClassificationContextWrapper&) = delete;

    ClassificationContextWrapper(ClassificationContextWrapper&& other) noexcept
        : ctx_(other.ctx_) {
        other.ctx_ = nullptr;
    }

    ClassificationContextWrapper& operator=(ClassificationContextWrapper&& other) noexcept {
        if (this != &other) {
            if (ctx_) {
                LLKA_destroyClassificationContext(ctx_);
            }
            ctx_ = other.ctx_;
            other.ctx_ = nullptr;
        }
        return *this;
    }

    LLKA_ClassificationContext* get() const { return ctx_; }
};

// Python bindings module
PYBIND11_MODULE(_pyllka_core, m) {
    m.doc() = "pyllka - Python bindings for libLLKA";

    // Register exception
    py::register_exception<LLKAError>(m, "LLKAError");

    // Enums
    py::enum_<LLKA_RetCode>(m, "RetCode")
        .value("OK", LLKA_OK)
        .value("INVALID_ARGUMENT", LLKA_E_INVALID_ARGUMENT)
        .value("MISMATCHING_SIZES", LLKA_E_MISMATCHING_SIZES)
        .value("NOT_IMPLEMENTED", LLKA_E_NOT_IMPLEMENTED)
        .value("MULTIPLE_ALT_IDS", LLKA_E_MULTIPLE_ALT_IDS)
        .value("MISSING_ATOMS", LLKA_E_MISSING_ATOMS)
        .value("MISMATCHING_DATA", LLKA_E_MISMATCHING_DATA)
        .value("BAD_CLASSIFICATION_CLUSTERS", LLKA_E_BAD_CLASSIFICATION_CLUSTERS)
        .value("BAD_GOLDEN_STEPS", LLKA_E_BAD_GOLDEN_STEPS)
        .value("BAD_CONFALS", LLKA_E_BAD_CONFALS)
        .value("BAD_AVERAGE_NU_ANGLES", LLKA_E_BAD_AVERAGE_NU_ANGLES)
        .value("BAD_CLASSIFICATION_LIMITS", LLKA_E_BAD_CLASSIFICATION_LIMITS)
        .value("NO_FILE", LLKA_E_NO_FILE)
        .value("CANNOT_READ_FILE", LLKA_E_CANNOT_READ_FILE)
        .value("BAD_DATA", LLKA_E_BAD_DATA)
        .value("NO_DATA", LLKA_E_NO_DATA)
        .value("NOTHING_TO_CLASSIFY", LLKA_E_NOTHING_TO_CLASSIFY)
        .export_values();

    py::enum_<LLKA_NtC>(m, "NtC")
        // AA classes
        .value("AA00", LLKA_AA00).value("AA01", LLKA_AA01).value("AA02", LLKA_AA02)
        .value("AA03", LLKA_AA03).value("AA04", LLKA_AA04).value("AA05", LLKA_AA05)
        .value("AA06", LLKA_AA06).value("AA07", LLKA_AA07).value("AA08", LLKA_AA08)
        .value("AA09", LLKA_AA09).value("AA10", LLKA_AA10).value("AA11", LLKA_AA11)
        .value("AA12", LLKA_AA12).value("AA13", LLKA_AA13)
        // AB classes
        .value("AB01", LLKA_AB01).value("AB02", LLKA_AB02).value("AB03", LLKA_AB03)
        .value("AB04", LLKA_AB04).value("AB05", LLKA_AB05)
        // BA classes
        .value("BA01", LLKA_BA01).value("BA05", LLKA_BA05).value("BA08", LLKA_BA08)
        .value("BA09", LLKA_BA09).value("BA10", LLKA_BA10).value("BA13", LLKA_BA13)
        .value("BA16", LLKA_BA16).value("BA17", LLKA_BA17)
        // BB classes
        .value("BB00", LLKA_BB00).value("BB01", LLKA_BB01).value("BB02", LLKA_BB02)
        .value("BB03", LLKA_BB03).value("BB04", LLKA_BB04).value("BB05", LLKA_BB05)
        .value("BB07", LLKA_BB07).value("BB08", LLKA_BB08).value("BB10", LLKA_BB10)
        .value("BB11", LLKA_BB11).value("BB12", LLKA_BB12).value("BB13", LLKA_BB13)
        .value("BB14", LLKA_BB14).value("BB15", LLKA_BB15).value("BB16", LLKA_BB16)
        .value("BB17", LLKA_BB17).value("BB20", LLKA_BB20)
        // IC classes
        .value("IC01", LLKA_IC01).value("IC02", LLKA_IC02).value("IC03", LLKA_IC03)
        .value("IC04", LLKA_IC04).value("IC05", LLKA_IC05).value("IC06", LLKA_IC06)
        .value("IC07", LLKA_IC07)
        // OP classes
        .value("OP01", LLKA_OP01).value("OP02", LLKA_OP02).value("OP03", LLKA_OP03)
        .value("OP04", LLKA_OP04).value("OP05", LLKA_OP05).value("OP06", LLKA_OP06)
        .value("OP07", LLKA_OP07).value("OP08", LLKA_OP08).value("OP09", LLKA_OP09)
        .value("OP10", LLKA_OP10).value("OP11", LLKA_OP11).value("OP12", LLKA_OP12)
        .value("OP13", LLKA_OP13).value("OP14", LLKA_OP14).value("OP15", LLKA_OP15)
        .value("OP16", LLKA_OP16).value("OP17", LLKA_OP17).value("OP18", LLKA_OP18)
        .value("OP19", LLKA_OP19).value("OP20", LLKA_OP20).value("OP21", LLKA_OP21)
        .value("OP22", LLKA_OP22).value("OP23", LLKA_OP23).value("OP24", LLKA_OP24)
        .value("OP25", LLKA_OP25).value("OP26", LLKA_OP26).value("OP27", LLKA_OP27)
        .value("OP28", LLKA_OP28).value("OP29", LLKA_OP29).value("OP30", LLKA_OP30)
        .value("OP31", LLKA_OP31).value("OPS1", LLKA_OPS1).value("OP1S", LLKA_OP1S)
        // Syn classes
        .value("AAS1", LLKA_AAS1).value("AB1S", LLKA_AB1S).value("AB2S", LLKA_AB2S)
        .value("BB1S", LLKA_BB1S).value("BB2S", LLKA_BB2S).value("BBS1", LLKA_BBS1)
        // ZZ classes
        .value("ZZ01", LLKA_ZZ01).value("ZZ02", LLKA_ZZ02)
        .value("ZZ1S", LLKA_ZZ1S).value("ZZ2S", LLKA_ZZ2S)
        .value("ZZS1", LLKA_ZZS1).value("ZZS2", LLKA_ZZS2)
        .value("INVALID_NTC", LLKA_INVALID_NTC)
        .export_values();

    py::enum_<LLKA_CANA>(m, "CANA")
        .value("AAA", LLKA_AAA).value("AAw", LLKA_AAw).value("AAu", LLKA_AAu)
        .value("AB", LLKA_AB)
        .value("BA", LLKA_BA)
        .value("BBB", LLKA_BBB).value("BBw", LLKA_BBw).value("B12", LLKA_B12)
        .value("BB2", LLKA_BB2).value("miB", LLKA_miB)
        .value("ICL", LLKA_ICL)
        .value("OPN", LLKA_OPN)
        .value("SYN", LLKA_SYN)
        .value("ZZZ", LLKA_ZZZ)
        .value("INVALID_CANA", LLKA_INVALID_CANA)
        .export_values();

    py::enum_<LLKA_SugarPucker>(m, "SugarPucker")
        .value("C3_ENDO", LLKA_C3_ENDO)
        .value("C4_EXO", LLKA_C4_EXO)
        .value("O4_ENDO", LLKA_O4_ENDO)
        .value("C1_EXO", LLKA_C1_EXO)
        .value("C2_ENDO", LLKA_C2_ENDO)
        .value("INVALID_SUGAR_PUCKER", LLKA_INVALID_SUGAR_PUCKER)
        .export_values();

    py::enum_<LLKA_SugarPuckerNameBrevity>(m, "SugarPuckerNameBrevity")
        .value("VERY_TERSE", LLKA_SPN_VERY_TERSE)
        .value("TERSE", LLKA_SPN_TERSE)
        .value("FANCY", LLKA_SPN_FANCY)
        .export_values();

    // Classification violation flags
    py::setattr(m, "CLASSIFICATION_OK", py::int_(static_cast<int>(LLKA_CLASSIFICATION_OK)));
    py::setattr(m, "CLASSIFICATION_E_SCORE_TOO_LOW", py::int_(static_cast<int>(LLKA_CLASSIFICATION_E_SCORE_TOO_LOW)));
    py::setattr(m, "CLASSIFICATION_E_NOT_ENOUGH_NEAREST_NEIGHBORS", py::int_(static_cast<int>(LLKA_CLASSIFICATION_E_NOT_ENOUGH_NEAREST_NEIGHBORS)));
    py::setattr(m, "CLASSIFICATION_E_AVERAGE_NEAREST_NEIGHBORS_TORSIONS_TOO_DIFFERENT", py::int_(static_cast<int>(LLKA_CLASSIFICATION_E_AVERAGE_NEAREST_NEIGHBORS_TORSIONS_TOO_DIFFERENT)));
    py::setattr(m, "CLASSIFICATION_E_NEAREST_NEIGHBOR_TORSIONS_TOO_DIFFERENT", py::int_(static_cast<int>(LLKA_CLASSIFICATION_E_NEAREST_NEIGHBOR_TORSIONS_TOO_DIFFERENT)));
    py::setattr(m, "CLASSIFICATION_E_CC_TOO_LOW", py::int_(static_cast<int>(LLKA_CLASSIFICATION_E_CC_TOO_LOW)));
    py::setattr(m, "CLASSIFICATION_E_CC_TOO_HIGH", py::int_(static_cast<int>(LLKA_CLASSIFICATION_E_CC_TOO_HIGH)));
    py::setattr(m, "CLASSIFICATION_E_NN_TOO_LOW", py::int_(static_cast<int>(LLKA_CLASSIFICATION_E_NN_TOO_LOW)));
    py::setattr(m, "CLASSIFICATION_E_NN_TOO_HIGH", py::int_(static_cast<int>(LLKA_CLASSIFICATION_E_NN_TOO_HIGH)));
    py::setattr(m, "CLASSIFICATION_E_MU_TOO_LOW", py::int_(static_cast<int>(LLKA_CLASSIFICATION_E_MU_TOO_LOW)));
    py::setattr(m, "CLASSIFICATION_E_MU_TOO_HIGH", py::int_(static_cast<int>(LLKA_CLASSIFICATION_E_MU_TOO_HIGH)));
    py::setattr(m, "CLASSIFICATION_E_TOTAL_DISTANCE_TOO_HIGH", py::int_(static_cast<int>(LLKA_CLASSIFICATION_E_TOTAL_DISTANCE_TOO_HIGH)));
    py::setattr(m, "CLASSIFICATION_E_FIRST_PSEUDOROTATION_TOO_DIFFERENT", py::int_(static_cast<int>(LLKA_CLASSIFICATION_E_FIRST_PSEUDOROTATION_TOO_DIFFERENT)));
    py::setattr(m, "CLASSIFICATION_E_SECOND_PSEUDOROTATION_TOO_DIFFERENT", py::int_(static_cast<int>(LLKA_CLASSIFICATION_E_SECOND_PSEUDOROTATION_TOO_DIFFERENT)));
    py::setattr(m, "CLASSIFICATION_E_BEST_CLUSTER_DOES_NOT_HAVE_ENOUGH_VOTES", py::int_(static_cast<int>(LLKA_CLASSIFICATION_E_BEST_CLUSTER_DOES_NOT_HAVE_ENOUGH_VOTES)));
    py::setattr(m, "CLASSIFICATION_E_DELTA_TORSION_ANGLE_REJECTED", py::int_(static_cast<int>(LLKA_CLASSIFICATION_E_DELTA_TORSION_ANGLE_REJECTED)));
    py::setattr(m, "CLASSIFICATION_E_WRONG_METRICS", py::int_(static_cast<int>(LLKA_CLASSIFICATION_E_WRONG_METRICS)));
    py::setattr(m, "CLASSIFICATION_E_UNASSIGNED_BUT_CLOSE_ENOUGH", py::int_(static_cast<int>(LLKA_CLASSIFICATION_E_UNASSIGNED_BUT_CLOSE_ENOUGH)));

    // Point
    py::class_<LLKA_Point>(m, "Point")
        .def(py::init<>())
        .def(py::init<double, double, double>())
        .def_readwrite("x", &LLKA_Point::x)
        .def_readwrite("y", &LLKA_Point::y)
        .def_readwrite("z", &LLKA_Point::z)
        .def("__repr__", [](const LLKA_Point& p) {
            return "Point(" + std::to_string(p.x) + ", " +
                   std::to_string(p.y) + ", " + std::to_string(p.z) + ")";
        });

    // Atom (read-only view)
    py::class_<LLKA_Atom>(m, "Atom")
        .def_property_readonly("type_symbol", [](const LLKA_Atom& a) {
            return std::string(a.type_symbol ? a.type_symbol : "");
        })
        .def_property_readonly("label_atom_id", [](const LLKA_Atom& a) {
            return std::string(a.label_atom_id ? a.label_atom_id : "");
        })
        .def_property_readonly("label_comp_id", [](const LLKA_Atom& a) {
            return std::string(a.label_comp_id ? a.label_comp_id : "");
        })
        .def_property_readonly("label_asym_id", [](const LLKA_Atom& a) {
            return std::string(a.label_asym_id ? a.label_asym_id : "");
        })
        .def_property_readonly("label_entity_id", [](const LLKA_Atom& a) {
            return std::string(a.label_entity_id ? a.label_entity_id : "");
        })
        .def_property_readonly("auth_atom_id", [](const LLKA_Atom& a) {
            return std::string(a.auth_atom_id ? a.auth_atom_id : "");
        })
        .def_property_readonly("auth_comp_id", [](const LLKA_Atom& a) {
            return std::string(a.auth_comp_id ? a.auth_comp_id : "");
        })
        .def_property_readonly("auth_asym_id", [](const LLKA_Atom& a) {
            return std::string(a.auth_asym_id ? a.auth_asym_id : "");
        })
        .def_property_readonly("pdbx_PDB_ins_code", [](const LLKA_Atom& a) {
            return std::string(a.pdbx_PDB_ins_code ? a.pdbx_PDB_ins_code : "");
        })
        .def_readonly("coords", &LLKA_Atom::coords)
        .def_readonly("id", &LLKA_Atom::id)
        .def_readonly("label_seq_id", &LLKA_Atom::label_seq_id)
        .def_readonly("auth_seq_id", &LLKA_Atom::auth_seq_id)
        .def_readonly("pdbx_PDB_model_num", &LLKA_Atom::pdbx_PDB_model_num)
        .def_readonly("label_alt_id", &LLKA_Atom::label_alt_id)
        .def("__repr__", [](const LLKA_Atom& a) {
            return "Atom(" + std::string(a.label_atom_id ? a.label_atom_id : "?") +
                   ", " + std::string(a.label_comp_id ? a.label_comp_id : "?") +
                   ", " + std::to_string(a.label_seq_id) + ")";
        });

    // Structure
    py::class_<StructureWrapper>(m, "Structure")
        .def(py::init<>())
        .def_property_readonly("n_atoms", &StructureWrapper::n_atoms)
        .def("get_atom", &StructureWrapper::get_atom, py::return_value_policy::reference_internal)
        .def("get_coordinates", &StructureWrapper::get_coordinates)
        .def("update_coordinates", &StructureWrapper::update_coordinates,
             py::arg("coords"),
             "Update structure coordinates from numpy array (n_atoms, 3)")
        .def("__len__", &StructureWrapper::n_atoms)
        .def("__getitem__", &StructureWrapper::get_atom, py::return_value_policy::reference_internal)
        .def("__repr__", [](const StructureWrapper& s) {
            return "Structure(n_atoms=" + std::to_string(s.n_atoms()) + ")";
        });

    // ImportedStructure (holds structure + CIF data)
    py::class_<ImportedStructureWrapper>(m, "ImportedStructure")
        .def(py::init<>())
        .def("get_structure", &ImportedStructureWrapper::get_structure,
             "Get the structure (creates a copy)")
        .def("get_entry_id", &ImportedStructureWrapper::get_entry_id,
             "Get entry ID from CIF")
        .def("__repr__", [](const ImportedStructureWrapper& is) {
            return "ImportedStructure(entry_id='" + is.get_entry_id() + "')";
        });

    // ClassificationContext
    py::class_<ClassificationContextWrapper>(m, "ClassificationContext")
        .def(py::init<>())
        .def("__repr__", [](const ClassificationContextWrapper&) {
            return "ClassificationContext()";
        });

    // StepMetrics
    // NuAngles (ribose torsion angles)
    py::class_<LLKA_NuAngles>(m, "NuAngles")
        .def_readonly("nu_0", &LLKA_NuAngles::nu_0)
        .def_readonly("nu_1", &LLKA_NuAngles::nu_1)
        .def_readonly("nu_2", &LLKA_NuAngles::nu_2)
        .def_readonly("nu_3", &LLKA_NuAngles::nu_3)
        .def_readonly("nu_4", &LLKA_NuAngles::nu_4);

    py::class_<LLKA_StepMetrics>(m, "StepMetrics")
        .def_readonly("delta_1", &LLKA_StepMetrics::delta_1)
        .def_readonly("epsilon_1", &LLKA_StepMetrics::epsilon_1)
        .def_readonly("zeta_1", &LLKA_StepMetrics::zeta_1)
        .def_readonly("alpha_2", &LLKA_StepMetrics::alpha_2)
        .def_readonly("beta_2", &LLKA_StepMetrics::beta_2)
        .def_readonly("gamma_2", &LLKA_StepMetrics::gamma_2)
        .def_readonly("delta_2", &LLKA_StepMetrics::delta_2)
        .def_readonly("chi_1", &LLKA_StepMetrics::chi_1)
        .def_readonly("chi_2", &LLKA_StepMetrics::chi_2)
        .def_readonly("CC", &LLKA_StepMetrics::CC)
        .def_readonly("NN", &LLKA_StepMetrics::NN)
        .def_readonly("mu", &LLKA_StepMetrics::mu);

    // ConfalScore
    py::class_<LLKA_ConfalScore>(m, "ConfalScore")
        .def_readonly("delta_1", &LLKA_ConfalScore::delta_1)
        .def_readonly("epsilon_1", &LLKA_ConfalScore::epsilon_1)
        .def_readonly("zeta_1", &LLKA_ConfalScore::zeta_1)
        .def_readonly("alpha_2", &LLKA_ConfalScore::alpha_2)
        .def_readonly("beta_2", &LLKA_ConfalScore::beta_2)
        .def_readonly("gamma_2", &LLKA_ConfalScore::gamma_2)
        .def_readonly("delta_2", &LLKA_ConfalScore::delta_2)
        .def_readonly("chi_1", &LLKA_ConfalScore::chi_1)
        .def_readonly("chi_2", &LLKA_ConfalScore::chi_2)
        .def_readonly("CC", &LLKA_ConfalScore::CC)
        .def_readonly("NN", &LLKA_ConfalScore::NN)
        .def_readonly("mu", &LLKA_ConfalScore::mu)
        .def_readonly("total", &LLKA_ConfalScore::total);

    // AverageConfal - average confal score and percentile
    py::class_<LLKA_AverageConfal>(m, "AverageConfal")
        .def_readonly("score", &LLKA_AverageConfal::score)
        .def_readonly("percentile", &LLKA_AverageConfal::percentile);

    // ClassifiedStep result
    py::class_<LLKA_ClassifiedStep>(m, "ClassifiedStep")
        .def_readonly("assigned_ntc", &LLKA_ClassifiedStep::assignedNtC)
        .def_readonly("assigned_cana", &LLKA_ClassifiedStep::assignedCANA)
        .def_readonly("closest_ntc", &LLKA_ClassifiedStep::closestNtC)
        .def_readonly("closest_cana", &LLKA_ClassifiedStep::closestCANA)
        .def_readonly("confal_score", &LLKA_ClassifiedStep::confalScore)
        .def_readonly("euclidean_distance_ntc_ideal", &LLKA_ClassifiedStep::euclideanDistanceNtCIdeal)
        .def_readonly("metrics", &LLKA_ClassifiedStep::metrics)
        .def_readonly("differences_from_ntc_averages", &LLKA_ClassifiedStep::differencesFromNtCAverages)
        .def_readonly("nu_angles_1", &LLKA_ClassifiedStep::nuAngles_1)
        .def_readonly("nu_angles_2", &LLKA_ClassifiedStep::nuAngles_2)
        .def_readonly("ribose_pseudorotation_1", &LLKA_ClassifiedStep::ribosePseudorotation_1)
        .def_readonly("ribose_pseudorotation_2", &LLKA_ClassifiedStep::ribosePseudorotation_2)
        .def_readonly("tau_1", &LLKA_ClassifiedStep::tau_1)
        .def_readonly("tau_2", &LLKA_ClassifiedStep::tau_2)
        .def_readonly("sugar_pucker_1", &LLKA_ClassifiedStep::sugarPucker_1)
        .def_readonly("sugar_pucker_2", &LLKA_ClassifiedStep::sugarPucker_2)
        .def_readonly("nu_angle_differences_1", &LLKA_ClassifiedStep::nuAngleDifferences_1)
        .def_readonly("nu_angle_differences_2", &LLKA_ClassifiedStep::nuAngleDifferences_2)
        .def_readonly("rmsd_to_closest_ntc", &LLKA_ClassifiedStep::rmsdToClosestNtC)
        .def_property_readonly("closest_golden_step", [](const LLKA_ClassifiedStep& cs) {
            return cs.closestGoldenStep ? std::string(cs.closestGoldenStep) : std::string("");
        })
        .def_readonly("violations", &LLKA_ClassifiedStep::violations)
        .def_readonly("violating_torsions_average", &LLKA_ClassifiedStep::violatingTorsionsAverage)
        .def_readonly("violating_torsions_nearest", &LLKA_ClassifiedStep::violatingTorsionsNearest)
        .def("__repr__", [](const LLKA_ClassifiedStep& cs) {
            return "ClassifiedStep(ntc=" + std::string(LLKA_NtCToName(cs.assignedNtC)) +
                   ", cana=" + std::string(LLKA_CANAToName(cs.assignedCANA)) +
                   ", rmsd=" + std::to_string(cs.rmsdToClosestNtC) + ")";
        });

    // I/O functions
    m.def("load_structure_from_file", [](const std::string& path) {
        LLKA_ImportedStructure imported;
        char* error = nullptr;

        LLKA_RetCode ret = LLKA_cifFileToStructure(
            path.c_str(), &imported, &error, LLKA_FALSE
        );

        if (ret != LLKA_OK) {
            std::string err_msg = error ? std::string(error) : "Unknown error";
            if (error) LLKA_destroyString(error);
            check_retcode(ret, ("Loading structure from " + path).c_str());
        }

        // Move structure to wrapper (transfer ownership)
        StructureWrapper wrapper(imported.structure);

        // Don't need to destroy imported since we moved the structure
        return wrapper;
    }, py::arg("path"), "Load structure from mmCIF file");

    // Load structure WITH CIF data (for writing back with annotations)
    m.def("load_structure_with_cif", [](const std::string& path) {
        LLKA_ImportedStructure imported;
        char* error = nullptr;

        LLKA_RetCode ret = LLKA_cifFileToStructure(
            path.c_str(), &imported, &error, LLKA_MINICIF_GET_CIFDATA
        );

        if (ret != LLKA_OK) {
            std::string err_msg = error ? std::string(error) : "Unknown error";
            if (error) LLKA_destroyString(error);
            check_retcode(ret, ("Loading structure with CIF data from " + path).c_str());
        }

        return ImportedStructureWrapper(imported);
    }, py::arg("path"), "Load structure with CIF data from mmCIF file (keeps original CIF for writing back)");

    // Write CIF data to string
    m.def("cif_to_string", [](const ImportedStructureWrapper& imported, bool pretty) {
        char* cifString = nullptr;
        LLKA_RetCode ret = LLKA_cifDataToString(imported.get_cif_data(), pretty ? LLKA_TRUE : LLKA_FALSE, &cifString);

        if (ret != LLKA_OK) {
            check_retcode(ret, "Converting CIF data to string");
        }

        std::string result(cifString);
        LLKA_destroyString(cifString);
        return result;
    }, py::arg("imported_structure"), py::arg("pretty") = true, "Convert CIF data to string");

    // Resource loading
    m.def("load_resource_file", [](const std::string& path) {
        LLKA_Resource resource;
        LLKA_RetCode ret = LLKA_loadResourceFile(path.c_str(), &resource);
        check_retcode(ret, ("Loading resource from " + path).c_str());
        return resource;
    }, py::arg("path"), "Load resource (CSV) file");

    // Utility functions
    m.def("ntc_to_name", [](LLKA_NtC ntc) {
        return std::string(LLKA_NtCToName(ntc));
    }, py::arg("ntc"), "Convert NtC enum to string name");

    m.def("name_to_ntc", [](const std::string& name) {
        return LLKA_nameToNtC(name.c_str());
    }, py::arg("name"), "Convert name string to NtC enum");

    m.def("cana_to_name", [](LLKA_CANA cana) {
        return std::string(LLKA_CANAToName(cana));
    }, py::arg("cana"), "Convert CANA enum to string name");

    m.def("sugar_pucker_to_name", [](LLKA_SugarPucker pucker, LLKA_SugarPuckerNameBrevity brevity = LLKA_SPN_TERSE) {
        return std::string(LLKA_sugarPuckerToName(pucker, brevity));
    }, py::arg("pucker"), py::arg("brevity") = LLKA_SPN_TERSE, "Convert SugarPucker enum to string name");

    m.def("classification_violation_to_name", [](int32_t violation) {
        return std::string(LLKA_classificationViolationToName(violation));
    }, py::arg("violation"), "Convert classification violation flag to string name");

    // Resource loading with type specification
    m.def("load_resource_file_typed", [](const std::string& path, LLKA_ResourceType type) {
        LLKA_Resource resource;
        resource.type = type;
        LLKA_RetCode ret = LLKA_loadResourceFile(path.c_str(), &resource);
        check_retcode(ret, ("Loading resource from " + path).c_str());

        // Return as Python dict with type and count
        py::dict result;
        result["type"] = static_cast<int>(type);
        result["count"] = resource.count;
        result["data"] = py::cast(&resource, py::return_value_policy::take_ownership);
        return result;
    }, py::arg("path"), py::arg("type"), "Load resource file with specific type");

    // Classification context creation
    m.def("create_classification_context_impl",
        [](const std::string& clusters_path,
           const std::string& golden_steps_path,
           const std::string& confals_path,
           const std::string& nu_angles_path,
           const std::string& confal_percentiles_path) -> ClassificationContextWrapper {

        // Load all resources
        LLKA_Resource clusters, goldenSteps, confals, nuAngles, confalPercentiles;

        clusters.type = LLKA_RES_CLUSTERS;
        LLKA_RetCode ret = LLKA_loadResourceFile(clusters_path.c_str(), &clusters);
        check_retcode(ret, "Loading clusters");

        goldenSteps.type = LLKA_RES_GOLDEN_STEPS;
        ret = LLKA_loadResourceFile(golden_steps_path.c_str(), &goldenSteps);
        if (ret != LLKA_OK) {
            LLKA_destroyResource(&clusters);
            check_retcode(ret, "Loading golden steps");
        }

        confals.type = LLKA_RES_CONFALS;
        ret = LLKA_loadResourceFile(confals_path.c_str(), &confals);
        if (ret != LLKA_OK) {
            LLKA_destroyResource(&clusters);
            LLKA_destroyResource(&goldenSteps);
            check_retcode(ret, "Loading confals");
        }

        nuAngles.type = LLKA_RES_AVERAGE_NU_ANGLES;
        ret = LLKA_loadResourceFile(nu_angles_path.c_str(), &nuAngles);
        if (ret != LLKA_OK) {
            LLKA_destroyResource(&clusters);
            LLKA_destroyResource(&goldenSteps);
            LLKA_destroyResource(&confals);
            check_retcode(ret, "Loading nu angles");
        }

        confalPercentiles.type = LLKA_RES_CONFAL_PERCENTILES;
        ret = LLKA_loadResourceFile(confal_percentiles_path.c_str(), &confalPercentiles);
        if (ret != LLKA_OK) {
            LLKA_destroyResource(&clusters);
            LLKA_destroyResource(&goldenSteps);
            LLKA_destroyResource(&confals);
            LLKA_destroyResource(&nuAngles);
            check_retcode(ret, "Loading confal percentiles");
        }

        // Set up classification limits (default values matching C++ examples)
        LLKA_ClassificationLimits limits;
        limits.minimumNearestNeighbors = 7;
        limits.numberOfUsedNearestNeighbors = 11;
        limits.averageNeighborsTorsionCutoff = 0.48869219055841207;  // 28 degrees in radians
        limits.nearestNeighborTorsionsCutoff = 0.48869219055841207;  // 28 degrees in radians
        limits.totalDistanceCutoff = 1.0471975511965976;  // 60 degrees in radians
        limits.pseudorotationCutoff = 1.2566370614359172;  // 72 degrees in radians
        limits.minimumClusterVotes = 0.001111;

        double maxCloseEnoughRmsd = 0.5;  // Angstroms

        // Create classification context
        LLKA_ClassificationContext* ctx = nullptr;
        ret = LLKA_initializeClassificationContext(
            clusters.data.clusters, clusters.count,
            goldenSteps.data.goldenSteps, goldenSteps.count,
            confals.data.confals, confals.count,
            nuAngles.data.clusterNuAngles, nuAngles.count,
            confalPercentiles.data.confalPercentiles, confalPercentiles.count,
            &limits,
            maxCloseEnoughRmsd,
            &ctx
        );

        // Clean up resources (context makes its own copies)
        LLKA_destroyResource(&clusters);
        LLKA_destroyResource(&goldenSteps);
        LLKA_destroyResource(&confals);
        LLKA_destroyResource(&nuAngles);
        LLKA_destroyResource(&confalPercentiles);

        check_retcode(ret, "Creating classification context");

        return ClassificationContextWrapper(ctx);
    },
    py::arg("clusters_path"),
    py::arg("golden_steps_path"),
    py::arg("confals_path"),
    py::arg("nu_angles_path"),
    py::arg("confal_percentiles_path"),
    "Create classification context from CSV files");

    // Structure splitting
    m.def("split_to_dinucleotides_impl", [](const StructureWrapper& structure) {
        LLKA_Structures steps;
        LLKA_RetCode ret = LLKA_splitStructureToDinucleotideSteps(structure.get(), &steps);
        check_retcode(ret, "Splitting structure to dinucleotides");

        // Convert to Python list of StructureWrapper
        py::list result;
        for (size_t i = 0; i < steps.nStrus; ++i) {
            result.append(StructureWrapper(steps.strus[i]));
        }

        // Free the array (but not the structures, as they're now owned by wrappers)
        free(steps.strus);

        return result;
    }, py::arg("structure"), "Split structure into dinucleotide steps");

    // Single step classification
    m.def("classify_step_impl", [](const StructureWrapper& structure,
                                     const ClassificationContextWrapper& context) {
        LLKA_ClassifiedStep result;
        LLKA_RetCode ret = LLKA_classifyStep(structure.get(), context.get(), &result);
        check_retcode(ret, "Classifying step");
        return result;
    }, py::arg("structure"), py::arg("context"), "Classify a single dinucleotide step");

    // Multiple steps classification
    m.def("classify_steps_multiple_impl",
        [](const py::list& structures, const ClassificationContextWrapper& context) {
        // Convert Python list to LLKA_Structures
        LLKA_Structures strus;
        strus.nStrus = structures.size();
        strus.strus = (LLKA_Structure*)malloc(sizeof(LLKA_Structure) * strus.nStrus);

        for (size_t i = 0; i < strus.nStrus; ++i) {
            const StructureWrapper* wrapper = structures[i].cast<const StructureWrapper*>();
            strus.strus[i] = *wrapper->get();
        }

        // Classify
        LLKA_ClassifiedSteps classifiedSteps;
        LLKA_RetCode ret = LLKA_classifyStepsMultiple(&strus, context.get(), &classifiedSteps);

        free(strus.strus);  // Free our temporary array

        check_retcode(ret, "Classifying multiple steps");

        // Convert results to Python list
        py::list results;
        for (size_t i = 0; i < classifiedSteps.nAttemptedSteps; ++i) {
            const auto& attempted = classifiedSteps.attemptedSteps[i];
            if (attempted.status == LLKA_OK) {
                results.append(attempted.step);
            } else {
                // Return step with violations
                results.append(attempted.step);
            }
        }

        LLKA_destroyClassifiedSteps(&classifiedSteps);

        return results;
    }, py::arg("structures"), py::arg("context"), "Classify multiple dinucleotide steps");

    // Average confal score and percentile
    m.def("average_confal_attempted", [](const py::list& classified_steps, const ClassificationContextWrapper& context) {
        // Convert Python list to LLKA_ClassifiedSteps
        LLKA_ClassifiedSteps steps;
        steps.nAttemptedSteps = classified_steps.size();
        steps.attemptedSteps = (LLKA_AttemptedClassifiedStep*)malloc(sizeof(LLKA_AttemptedClassifiedStep) * steps.nAttemptedSteps);

        for (size_t i = 0; i < steps.nAttemptedSteps; i++) {
            steps.attemptedSteps[i].step = classified_steps[i].cast<LLKA_ClassifiedStep>();
            steps.attemptedSteps[i].status = LLKA_OK;  // Assume OK status for steps that succeeded
        }

        LLKA_AverageConfal result = LLKA_averageConfalAttempted(&steps, context.get());

        free(steps.attemptedSteps);

        return result;
    }, py::arg("classified_steps"), py::arg("context"),
       "Calculate average confal score and percentile for a list of classified steps");

    // Version constants
    m.attr("NTC_VERSION") = LLKA_INTERNAL_NTC_VERSION;
    m.attr("CANA_VERSION") = LLKA_INTERNAL_CANA_VERSION;

    // Angle utility functions
    m.def("rad2deg", [](double radians) -> double {
        return LLKA_rad2deg(radians);
    }, py::arg("radians"), "Convert radians to degrees");

    m.def("full_angle_from_rad", [](double radians) -> double {
        return LLKA_fullAngleFromRad(radians);
    }, py::arg("radians"), "Normalize angle in radians to 0-2π range (adds 2π to negative angles)");

    // Measurement functions
    m.def("measure_distance",
        [](const LLKA_Atom* a, const LLKA_Atom* b) -> double {
            return LLKA_measureDistance(a, b);
        },
        py::arg("a"), py::arg("b"),
        "Measure distance between two atoms");

    m.def("measure_angle",
        [](const LLKA_Atom* a, const LLKA_Atom* b, const LLKA_Atom* c) -> double {
            return LLKA_measureAngle(a, b, c);
        },
        py::arg("a"), py::arg("b"), py::arg("c"),
        "Measure angle between three atoms (in radians)");

    m.def("measure_dihedral",
        [](const LLKA_Atom* a, const LLKA_Atom* b,
           const LLKA_Atom* c, const LLKA_Atom* d) -> double {
            return LLKA_measureDihedral(a, b, c, d);
        },
        py::arg("a"), py::arg("b"), py::arg("c"), py::arg("d"),
        "Measure dihedral angle between four atoms (in radians)");

    // Step metrics calculation
    m.def("calculate_step_metrics", [](const StructureWrapper& structure) {
        LLKA_StepMetrics metrics;
        LLKA_RetCode ret = LLKA_calculateStepMetrics(structure.get(), &metrics);
        check_retcode(ret, "Calculating step metrics");
        return metrics;
    }, py::arg("structure"), "Calculate NtC step metrics for a dinucleotide");

    // RMSD and superposition
    m.def("rmsd_structures", [](const StructureWrapper& a, const StructureWrapper& b) {
        double rmsd;
        LLKA_RetCode ret = LLKA_rmsdStructures(a.get(), b.get(), &rmsd);
        check_retcode(ret, "Calculating RMSD");
        return rmsd;
    }, py::arg("a"), py::arg("b"), "Calculate RMSD between two structures");

    m.def("superpose_structures", [](StructureWrapper& what, const StructureWrapper& onto) {
        double rmsd;
        LLKA_RetCode ret = LLKA_superposeStructures(what.get(), onto.get(), &rmsd);
        check_retcode(ret, "Superposing structures");
        return rmsd;
    }, py::arg("what"), py::arg("onto"),
       "Superpose 'what' onto 'onto' structure (modifies 'what' in place), returns RMSD");

    // Similarity and Connectivity
    py::class_<LLKA_Similarity>(m, "Similarity")
        .def(py::init<>())
        .def_readwrite("rmsd", &LLKA_Similarity::rmsd, "RMSD of the steps after Kabsch superposition")
        .def_readwrite("euclidean_distance", &LLKA_Similarity::euclideanDistance, "Sum of the euclidean distances between the torsions or distances")
        .def("__repr__", [](const LLKA_Similarity& s) {
            return "<Similarity rmsd=" + std::to_string(s.rmsd) +
                   " euclidean_distance=" + std::to_string(s.euclideanDistance) + ">";
        });

    py::class_<LLKA_Connectivity>(m, "Connectivity")
        .def(py::init<>())
        .def_readwrite("c5_prime_distance", &LLKA_Connectivity::C5PrimeDistance,
                      "Distance between C5' atoms")
        .def_readwrite("o3_prime_distance", &LLKA_Connectivity::O3PrimeDistance,
                      "Distance between O3' atoms")
        .def("__repr__", [](const LLKA_Connectivity& c) {
            return "<Connectivity c5_prime_distance=" + std::to_string(c.C5PrimeDistance) +
                   " o3_prime_distance=" + std::to_string(c.O3PrimeDistance) + ">";
        });

    // Similarity measurement functions
    m.def("measure_step_similarity_ntc",
        [](const StructureWrapper& step, LLKA_NtC ntc) {
            LLKA_Similarity result;
            LLKA_RetCode ret = LLKA_measureStepSimilarityNtC(step.get(), ntc, &result);
            check_retcode(ret, "Measuring step similarity");
            return result;
        },
        py::arg("step"), py::arg("ntc"),
        "Measure similarity between a dinucleotide and a reference NtC");

    m.def("measure_step_similarity_ntc_multiple",
        [](const StructureWrapper& step, const std::vector<LLKA_NtC>& ntcs) {
            // Create null-terminated array
            std::vector<LLKA_NtC> ntcs_array = ntcs;
            ntcs_array.push_back(LLKA_INVALID_NTC);

            // Pre-allocate results structure
            LLKA_Similarities results;
            results.nSimilars = ntcs.size();
            results.similars = (LLKA_Similarity*)malloc(sizeof(LLKA_Similarity) * results.nSimilars);
            if (!results.similars) {
                throw LLKAError("Failed to allocate memory for similarities");
            }

            LLKA_RetCode ret = LLKA_measureStepSimilarityNtCMultiple(step.get(), ntcs_array.data(), &results);
            if (ret != LLKA_OK) {
                free(results.similars);
                check_retcode(ret, "Measuring step similarity multiple");
            }

            // Convert to vector
            std::vector<LLKA_Similarity> result_vec(results.similars, results.similars + results.nSimilars);
            free(results.similars);

            return result_vec;
        },
        py::arg("step"), py::arg("ntcs"),
        "Measure similarity between a dinucleotide and multiple reference NtCs");

    m.def("measure_step_similarity_structure",
        [](const StructureWrapper& step, const StructureWrapper& reference) {
            LLKA_Similarity result;
            LLKA_RetCode ret = LLKA_measureStepSimilarityStructure(step.get(), reference.get(), &result);
            check_retcode(ret, "Measuring step similarity with structure");
            return result;
        },
        py::arg("step"), py::arg("reference"),
        "Measure similarity between two dinucleotides");

    // Connectivity measurement functions
    m.def("measure_step_connectivity_ntcs",
        [](const StructureWrapper& position_first, LLKA_NtC ntc_first,
           const StructureWrapper& position_second, LLKA_NtC ntc_second) {
            LLKA_Connectivity result;
            LLKA_RetCode ret = LLKA_measureStepConnectivityNtCs(
                position_first.get(), ntc_first,
                position_second.get(), ntc_second,
                &result);
            check_retcode(ret, "Measuring step connectivity");
            return result;
        },
        py::arg("position_first"), py::arg("ntc_first"),
        py::arg("position_second"), py::arg("ntc_second"),
        "Measure connectivity between two synthetic steps corresponding to NtCs");

    m.def("measure_step_connectivity_ntcs_multiple_first",
        [](const StructureWrapper& position_first, const std::vector<LLKA_NtC>& ntcs_first,
           const StructureWrapper& position_second, LLKA_NtC ntc_second) {
            // Create null-terminated array
            std::vector<LLKA_NtC> ntcs_array = ntcs_first;
            ntcs_array.push_back(LLKA_INVALID_NTC);

            // Pre-allocate results structure
            LLKA_Connectivities results;
            results.nConns = ntcs_first.size();
            results.conns = (LLKA_Connectivity*)malloc(sizeof(LLKA_Connectivity) * results.nConns);
            if (!results.conns) {
                throw LLKAError("Failed to allocate memory for connectivities");
            }

            LLKA_RetCode ret = LLKA_measureStepConnectivityNtCsMultipleFirst(
                position_first.get(), ntcs_array.data(),
                position_second.get(), ntc_second,
                &results);
            if (ret != LLKA_OK) {
                free(results.conns);
                check_retcode(ret, "Measuring step connectivity multiple first");
            }

            // Convert to vector
            std::vector<LLKA_Connectivity> result_vec(results.conns, results.conns + results.nConns);
            free(results.conns);

            return result_vec;
        },
        py::arg("position_first"), py::arg("ntcs_first"),
        py::arg("position_second"), py::arg("ntc_second"),
        "Measure connectivity with multiple NtCs for the first step");

    m.def("measure_step_connectivity_ntcs_multiple_second",
        [](const StructureWrapper& position_first, LLKA_NtC ntc_first,
           const StructureWrapper& position_second, const std::vector<LLKA_NtC>& ntcs_second) {
            // Create null-terminated array
            std::vector<LLKA_NtC> ntcs_array = ntcs_second;
            ntcs_array.push_back(LLKA_INVALID_NTC);

            // Pre-allocate results structure
            LLKA_Connectivities results;
            results.nConns = ntcs_second.size();
            results.conns = (LLKA_Connectivity*)malloc(sizeof(LLKA_Connectivity) * results.nConns);
            if (!results.conns) {
                throw LLKAError("Failed to allocate memory for connectivities");
            }

            LLKA_RetCode ret = LLKA_measureStepConnectivityNtCsMultipleSecond(
                position_first.get(), ntc_first,
                position_second.get(), ntcs_array.data(),
                &results);
            if (ret != LLKA_OK) {
                free(results.conns);
                check_retcode(ret, "Measuring step connectivity multiple second");
            }

            // Convert to vector
            std::vector<LLKA_Connectivity> result_vec(results.conns, results.conns + results.nConns);
            free(results.conns);

            return result_vec;
        },
        py::arg("position_first"), py::arg("ntc_first"),
        py::arg("position_second"), py::arg("ntcs_second"),
        "Measure connectivity with multiple NtCs for the second step");

    m.def("measure_step_connectivity_structures",
        [](const StructureWrapper& position_first, const StructureWrapper& dinu_first,
           const StructureWrapper& position_second, const StructureWrapper& dinu_second) {
            LLKA_Connectivity result;
            LLKA_RetCode ret = LLKA_measureStepConnectivityStructures(
                position_first.get(), dinu_first.get(),
                position_second.get(), dinu_second.get(),
                &result);
            check_retcode(ret, "Measuring step connectivity with structures");
            return result;
        },
        py::arg("position_first"), py::arg("dinu_first"),
        py::arg("position_second"), py::arg("dinu_second"),
        "Measure connectivity between two dinucleotides");
}
