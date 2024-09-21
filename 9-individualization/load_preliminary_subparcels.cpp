#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <memory>
#include "json.hpp" // Include the json.hpp header file

using json = nlohmann::json; // Define a shorter alias for nlohmann::json

// Define your data structures here
struct Vertex {
    int index;
    double x, y, z;
    std::string label_parcel;
    std::vector<int> neighbors_index;
    Vertex(int index, double x, double y, double z, std::string label_parcel) : index(index), x(x), y(y), z(z), label_parcel(label_parcel) {}
};

struct Triangle {
    int index;
    std::shared_ptr<Vertex> v1, v2, v3;
    int label_parcel;
    std::vector<int> label_subparcel;
    std::vector<int> fibers;
    std::unordered_set<std::shared_ptr<Triangle>> neighbors;
    Triangle(int index, std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2, std::shared_ptr<Vertex> v3, int label_parcel, const std::vector<int>& label_subparcel, const std::vector<int>& fibers) : index(index), v1(v1), v2(v2), v3(v3), label_parcel(label_parcel), label_subparcel(label_subparcel), fibers(fibers) {}
};

struct SubParcel {
    std::string label;
    std::string label_anatomic;
    std::vector<std::shared_ptr<Triangle>> triangles;
    std::vector<double> inter_points;
    std::vector<int> bundle_labels;
    SubParcel(std::string label, std::string label_anatomic, const std::vector<std::shared_ptr<Triangle>>& triangles, const std::vector<double>& inter_points) : label(label), label_anatomic(label_anatomic), triangles(triangles), inter_points(inter_points) {}
    void add_bundle_label(int label) {
        bundle_labels.push_back(label);
    }
};

struct AnatomicParcel {
    std::string label;
    std::unordered_map<std::string, std::shared_ptr<SubParcel>> sub_parcels;
    AnatomicParcel(std::string label) : label(label) {}
};

std::vector<AnatomicParcelData> parse_json(const std::string& json_string) {
    std::vector<AnatomicParcelData> data;

    // Parse JSON string
    json j = json::parse(json_string);

    // Extract data from JSON and populate AnatomicParcelData objects
    // Example code, replace with your actual JSON parsing logic
    for (const auto& ap : j["aparcels"]) {
    AnatomicParcelData ap_data;
    ap_data.label = ap["label"];

    // Process 'sub_parcels' array
    for (const auto& sp : ap["sub_parcels"]) {
        SubParcelData subparcel_data;
        subparcel_data.label = sp["label"];
        subparcel_data.label_anatomic = sp["label_anatomic"];

        // Process 'triangles' array
        for (const auto& tri : sp["triangles"]) {
            TriangleData triangle_data;
            triangle_data.index = tri["index"];
            // Extract and populate other fields of the triangle_data object

            subparcel_data.triangles.push_back(triangle_data);
        }

        // Process other fields of the SubParcelData object
        // subparcel_data.inter_points = ...
        // subparcel_data.bundle_labels = ...

        ap_data.sub_parcels.push_back(subparcel_data);
    }

    // Process other fields of the AnatomicParcelData object
    // ap_data.other_field = ...

    data.push_back(ap_data);
}

    return data;
}

std::vector<std::shared_ptr<AnatomicParcel>> load_preliminary_subparcels(const std::string& file_path) {
    std::ifstream file(file_path);
    if (!file.is_open()) {
        std::cerr << "Error opening file " << file_path << std::endl;
        return {};
    }

    // Read JSON string from file
    std::string json_string((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());

    // Parse JSON string and extract data
    std::vector<AnatomicParcelData> data = parse_json(json_string);

    std::vector<std::shared_ptr<AnatomicParcel>> aparcels;

    // Process data and create AnatomicParcel objects
    for (const auto& ap_data : data) {
        std::shared_ptr<AnatomicParcel> anatomic_parcel = std::make_shared<AnatomicParcel>(ap_data.label);
        
        for (const auto& subparcel_data : ap_data.sub_parcels) {
            std::vector<std::shared_ptr<Triangle>> triangles;
            for (const auto& triangle_data : subparcel_data.triangles) {
                // Create vertices
                std::shared_ptr<Vertex> v1 = std::make_shared<Vertex>(triangle_data.v1.index, triangle_data.v1.x, triangle_data.v1.y, triangle_data.v1.z, triangle_data.v1.label_parcel);
                std::shared_ptr<Vertex> v2 = std::make_shared<Vertex>(triangle_data.v2.index, triangle_data.v2.x, triangle_data.v2.y, triangle_data.v2.z, triangle_data.v2.label_parcel);
                std::shared_ptr<Vertex> v3 = std::make_shared<Vertex>(triangle_data.v3.index, triangle_data.v3.x, triangle_data.v3.y, triangle_data.v3.z, triangle_data.v3.label_parcel);

                // Create triangle
                std::shared_ptr<Triangle> triangle = std::make_shared<Triangle>(triangle_data.index, v1, v2, v3, triangle_data.label_parcel, triangle_data.label_subparcel, triangle_data.fibers);
                triangles.push_back(triangle);
            }
            // Create SubParcel
            std::shared_ptr<SubParcel> subparcel = std::make_shared<SubParcel>(subparcel_data.label, subparcel_data.label_anatomic, triangles, subparcel_data.inter_points);
            for (const auto& bundle_label : subparcel_data.bundle_labels) {
                subparcel->add_bundle_label(bundle_label);
            }

            anatomic_parcel->sub_parcels[subparcel->label] = subparcel;
        }
        aparcels.push_back(anatomic_parcel);
    }

    // Release memory after processing
    return aparcels;
}

int main() {
    // Example usage
    std::string file_path = "your_file.json";
    auto aparcels = load_preliminary_subparcels(file_path);
    // Use the loaded data
    return 0;
}