#include "ThreadsIO.hpp"

#include <hdf5.h>
#include <chrono>
#include <cstring>
#include <ctime>
#include <stdexcept>
#include <vector>

// ── helpers ──────────────────────────────────────────────────────────────────

namespace {

struct H5Handle {
    hid_t id;
    void (*closer)(hid_t);
    H5Handle(hid_t id, void (*closer)(hid_t)) : id(id), closer(closer) {
        if (id < 0) throw std::runtime_error("HDF5 handle creation failed");
    }
    ~H5Handle() { if (id >= 0) closer(id); }
    operator hid_t() const { return id; }
    H5Handle(const H5Handle&) = delete;
    H5Handle& operator=(const H5Handle&) = delete;
};

void close_file(hid_t h)    { H5Fclose(h); }
void close_space(hid_t h)   { H5Sclose(h); }
void close_dset(hid_t h)    { H5Dclose(h); }
void close_plist(hid_t h)   { H5Pclose(h); }
void close_type(hid_t h)    { H5Tclose(h); }
void close_attr(hid_t h)    { H5Aclose(h); }

hid_t make_dcpl(int rank, const hsize_t* dims) {
    // Cannot chunk/compress empty datasets
    for (int i = 0; i < rank; i++)
        if (dims[i] == 0) return H5P_DEFAULT;
    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl, rank, dims);
    H5Pset_deflate(dcpl, 9);
    return dcpl;
}

void write_i64_1d(hid_t file, const char* name, const std::vector<int64_t>& data) {
    hsize_t n = data.size();
    H5Handle space(H5Screate_simple(1, &n, nullptr), close_space);
    H5Handle dcpl(make_dcpl(1, &n), close_plist);
    H5Handle ds(H5Dcreate2(file, name, H5T_STD_I64LE, space, H5P_DEFAULT, dcpl, H5P_DEFAULT), close_dset);
    H5Dwrite(ds, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
}

void write_f64_1d(hid_t file, const char* name, const std::vector<double>& data) {
    hsize_t n = data.size();
    H5Handle space(H5Screate_simple(1, &n, nullptr), close_space);
    H5Handle dcpl(make_dcpl(1, &n), close_plist);
    H5Handle ds(H5Dcreate2(file, name, H5T_IEEE_F64LE, space, H5P_DEFAULT, dcpl, H5P_DEFAULT), close_dset);
    H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
}

void write_i64_2d(hid_t file, const char* name, const int64_t* data, hsize_t rows, hsize_t cols) {
    hsize_t dims[2] = {rows, cols};
    H5Handle space(H5Screate_simple(2, dims, nullptr), close_space);
    H5Handle dcpl(make_dcpl(2, dims), close_plist);
    H5Handle ds(H5Dcreate2(file, name, H5T_STD_I64LE, space, H5P_DEFAULT, dcpl, H5P_DEFAULT), close_dset);
    H5Dwrite(ds, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
}

void write_string_2d(hid_t file, const char* name, const std::vector<std::vector<std::string>>& cols, hsize_t rows) {
    hsize_t ncols = cols.size();
    hsize_t dims[2] = {rows, ncols};
    H5Handle space(H5Screate_simple(2, dims, nullptr), close_space);
    H5Handle dcpl(make_dcpl(2, dims), close_plist);
    H5Handle strtype(H5Tcopy(H5T_C_S1), close_type);
    H5Tset_size(strtype, H5T_VARIABLE);
    H5Handle ds(H5Dcreate2(file, name, strtype, space, H5P_DEFAULT, dcpl, H5P_DEFAULT), close_dset);

    // Build flat row-major array of const char*
    std::vector<const char*> ptrs(rows * ncols);
    for (hsize_t r = 0; r < rows; r++)
        for (hsize_t c = 0; c < ncols; c++)
            ptrs[r * ncols + c] = cols[c][r].c_str();
    H5Dwrite(ds, strtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, ptrs.data());
}

void write_string_1d(hid_t file, const char* name, const std::vector<std::string>& data) {
    hsize_t n = data.size();
    H5Handle space(H5Screate_simple(1, &n, nullptr), close_space);
    H5Handle dcpl(make_dcpl(1, &n), close_plist);
    H5Handle strtype(H5Tcopy(H5T_C_S1), close_type);
    H5Tset_size(strtype, H5T_VARIABLE);
    H5Handle ds(H5Dcreate2(file, name, strtype, space, H5P_DEFAULT, dcpl, H5P_DEFAULT), close_dset);
    std::vector<const char*> ptrs(n);
    for (hsize_t i = 0; i < n; i++) ptrs[i] = data[i].c_str();
    H5Dwrite(ds, strtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, ptrs.data());
}

std::vector<int64_t> read_i64(hid_t file, const char* name) {
    H5Handle ds(H5Dopen2(file, name, H5P_DEFAULT), close_dset);
    H5Handle space(H5Dget_space(ds), close_space);
    hsize_t dims[2];
    int ndims = H5Sget_simple_extent_dims(space, dims, nullptr);
    hsize_t total = dims[0];
    if (ndims == 2) total *= dims[1];
    std::vector<int64_t> data(total);
    H5Dread(ds, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
    return data;
}

std::vector<double> read_f64(hid_t file, const char* name) {
    H5Handle ds(H5Dopen2(file, name, H5P_DEFAULT), close_dset);
    H5Handle space(H5Dget_space(ds), close_space);
    hsize_t dims[2];
    int ndims = H5Sget_simple_extent_dims(space, dims, nullptr);
    hsize_t total = dims[0];
    if (ndims == 2) total *= dims[1];
    std::vector<double> data(total);
    H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
    return data;
}

std::vector<std::string> read_varlen_strings(hid_t file, const char* name) {
    H5Handle ds(H5Dopen2(file, name, H5P_DEFAULT), close_dset);
    H5Handle space(H5Dget_space(ds), close_space);
    hsize_t dims[2];
    int ndims = H5Sget_simple_extent_dims(space, dims, nullptr);
    hsize_t total = dims[0];
    if (ndims == 2) total *= dims[1];

    H5Handle memtype(H5Tcopy(H5T_C_S1), close_type);
    H5Tset_size(memtype, H5T_VARIABLE);

    std::vector<char*> raw(total);
    H5Dread(ds, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, raw.data());

    std::vector<std::string> result(total);
    for (hsize_t i = 0; i < total; i++) {
        result[i] = raw[i] ? raw[i] : "";
    }

    // Reclaim memory allocated by HDF5
    H5Handle filetype(H5Dget_type(ds), close_type);
    H5Dvlen_reclaim(memtype, space, H5P_DEFAULT, raw.data());
    return result;
}

hsize_t get_dim0(hid_t file, const char* name) {
    H5Handle ds(H5Dopen2(file, name, H5P_DEFAULT), close_dset);
    H5Handle space(H5Dget_space(ds), close_space);
    hsize_t dims[2];
    H5Sget_simple_extent_dims(space, dims, nullptr);
    return dims[0];
}

std::string iso_now() {
    auto now = std::chrono::system_clock::now();
    auto t = std::chrono::system_clock::to_time_t(now);
    char buf[64];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S", std::localtime(&t));
    return buf;
}

} // anonymous namespace

// ── public API ───────────────────────────────────────────────────────────────

void serialize_threads(
    const std::string& filename,
    ThreadingInstructions& instructions,
    const std::vector<std::vector<std::string>>& metadata_cols,
    const std::vector<double>& allele_ages,
    const std::vector<std::string>& sample_names)
{
    int N = instructions.num_samples;
    int M = instructions.num_sites;

    auto all_starts = instructions.all_starts();
    auto all_targets = instructions.all_targets();
    auto all_tmrcas = instructions.all_tmrcas();
    auto all_mismatches = instructions.all_mismatches();

    // Build samples array (N x 3) and flatten thread/mismatch data
    std::vector<int64_t> samples_flat(N * 3);
    std::vector<int64_t> flat_targets, flat_starts, flat_mismatches;
    std::vector<double> flat_tmrcas;

    int64_t toff = 0, moff = 0;
    for (int i = 0; i < N; i++) {
        samples_flat[i * 3 + 0] = i;
        samples_flat[i * 3 + 1] = toff;
        samples_flat[i * 3 + 2] = moff;
        for (auto v : all_targets[i])    flat_targets.push_back(v);
        for (auto v : all_starts[i])     flat_starts.push_back(v);
        for (auto v : all_tmrcas[i])     flat_tmrcas.push_back(v);
        for (auto v : all_mismatches[i]) flat_mismatches.push_back(v);
        toff += all_starts[i].size();
        moff += all_mismatches[i].size();
    }

    // Build thread_targets 2D array (S x 2): [target, start]
    size_t S = flat_targets.size();
    std::vector<int64_t> targets_2d(S * 2);
    for (size_t i = 0; i < S; i++) {
        targets_2d[i * 2 + 0] = flat_targets[i];
        targets_2d[i * 2 + 1] = flat_starts[i];
    }

    std::vector<double> arg_range = {
        static_cast<double>(instructions.start),
        static_cast<double>(instructions.end)
    };

    std::vector<int64_t> positions_i64(instructions.positions.begin(), instructions.positions.end());
    std::vector<int64_t> mismatches_i64(flat_mismatches.begin(), flat_mismatches.end());

    // Write file
    H5Handle file(H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT), close_file);

    // datetime attribute
    std::string dt = iso_now();
    {
        H5Handle strtype(H5Tcopy(H5T_C_S1), close_type);
        H5Tset_size(strtype, H5T_VARIABLE);
        H5Handle aspace(H5Screate(H5S_SCALAR), close_space);
        H5Handle attr(H5Acreate2(file, "datetime_created", strtype, aspace, H5P_DEFAULT, H5P_DEFAULT), close_attr);
        const char* dtp = dt.c_str();
        H5Awrite(attr, strtype, &dtp);
    }

    write_i64_2d(file, "samples", samples_flat.data(), N, 3);
    write_i64_1d(file, "positions", positions_i64);
    write_i64_2d(file, "thread_targets", targets_2d.data(), S, 2);
    write_f64_1d(file, "thread_ages", flat_tmrcas);
    write_i64_1d(file, "het_sites", mismatches_i64);
    write_f64_1d(file, "arg_range", arg_range);

    // Optional datasets
    if (!metadata_cols.empty() && M > 0) {
        write_string_2d(file, "variant_metadata", metadata_cols, M);
    }
    if (!allele_ages.empty()) {
        write_f64_1d(file, "allele_ages", allele_ages);
    }
    if (!sample_names.empty()) {
        write_string_1d(file, "sample_names", sample_names);
    }
}


ThreadingInstructions deserialize_threads(const std::string& filename) {
    H5Handle file(H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), close_file);

    auto samples_flat = read_i64(file, "samples");
    auto positions_i64 = read_i64(file, "positions");
    auto targets_flat = read_i64(file, "thread_targets");
    auto tmrcas_flat = read_f64(file, "thread_ages");
    auto mismatches_flat = read_i64(file, "het_sites");

    hsize_t N = get_dim0(file, "samples");

    // Read arg_range (may not exist in very old files)
    double region_start, region_end;
    if (H5Lexists(file, "arg_range", H5P_DEFAULT) > 0) {
        auto ar = read_f64(file, "arg_range");
        region_start = ar[0];
        region_end = ar[1];
    } else {
        region_start = 0;
        region_end = 0;
    }

    // Extract per-sample thread offsets and het offsets from samples array
    std::vector<int64_t> thread_starts(N), het_starts(N);
    for (hsize_t i = 0; i < N; i++) {
        thread_starts[i] = samples_flat[i * 3 + 1];
        het_starts[i] = samples_flat[i * 3 + 2];
    }

    hsize_t total_threads = targets_flat.size() / 2;
    hsize_t total_hets = mismatches_flat.size();

    // Split flat arrays into per-sample vectors
    std::vector<std::vector<int>> starts(N), targets(N), mismatches(N);
    std::vector<std::vector<double>> tmrcas(N);
    std::vector<int> positions(positions_i64.begin(), positions_i64.end());

    for (hsize_t i = 0; i < N; i++) {
        int64_t t_start = thread_starts[i];
        int64_t t_end = (i + 1 < N) ? thread_starts[i + 1] : static_cast<int64_t>(total_threads);
        int64_t h_start = het_starts[i];
        int64_t h_end = (i + 1 < N) ? het_starts[i + 1] : static_cast<int64_t>(total_hets);

        for (int64_t j = t_start; j < t_end; j++) {
            targets[i].push_back(static_cast<int>(targets_flat[j * 2]));
            starts[i].push_back(static_cast<int>(targets_flat[j * 2 + 1]));
            tmrcas[i].push_back(tmrcas_flat[j]);
        }
        for (int64_t j = h_start; j < h_end; j++) {
            mismatches[i].push_back(static_cast<int>(mismatches_flat[j]));
        }
    }

    return ThreadingInstructions(
        starts, tmrcas, targets, mismatches,
        positions,
        static_cast<int>(region_start),
        static_cast<int>(region_end)
    );
}


std::vector<std::vector<std::string>> read_threads_metadata(const std::string& filename) {
    H5Handle file(H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), close_file);
    if (H5Lexists(file, "variant_metadata", H5P_DEFAULT) <= 0) {
        throw std::runtime_error("No variant_metadata in " + filename);
    }

    hsize_t rows = get_dim0(file, "variant_metadata");
    auto flat = read_varlen_strings(file, "variant_metadata");

    // Unflatten: flat is row-major (rows x 7) → 7 columns
    const int ncols = 7;
    std::vector<std::vector<std::string>> cols(ncols);
    for (int c = 0; c < ncols; c++) {
        cols[c].resize(rows);
        for (hsize_t r = 0; r < rows; r++) {
            cols[c][r] = flat[r * ncols + c];
        }
    }
    return cols;
}


std::vector<std::string> read_threads_sample_names(const std::string& filename) {
    H5Handle file(H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), close_file);
    if (H5Lexists(file, "sample_names", H5P_DEFAULT) <= 0) {
        throw std::runtime_error("No sample_names in " + filename);
    }
    return read_varlen_strings(file, "sample_names");
}
