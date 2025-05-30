#include "VTKIOFunc.h"

namespace VTKIOFunc {
	void OutputCircleDisksAsVTP(const Array<Vector<real, 3>>& point_array, const Array<Vector<real, 3>>& normal_array, const Array<real>& area_array, const Array<real>& height_array, const fs::path& path) {
		Info("Outputting Lagrangian particles as disk glyphs to VTU file: {}", path.string());

		const int N = point_array.size();
		Assert(N == normal_array.size() && N == area_array.size(), "Points, normals, and areas must have the same size.");

		// Points
		vtkNew<vtkPoints> points;
		points->SetNumberOfPoints(N);
		for (int i = 0; i < N; ++i) {
			const auto& x = point_array[i];
			points->SetPoint(i, x[0], x[1], x[2]);
		}

		// Normals
		vtkNew<vtkFloatArray> normals;
		normals->SetNumberOfComponents(3);
		normals->SetNumberOfTuples(N);
		normals->SetName("Normals");
		for (int i = 0; i < N; ++i) {
			const auto& n = normal_array[i];
			normals->SetTuple3(i, n[0], n[1], n[2]);
		}

		// Scaling (radius based on area A, disk area = πr² => r = sqrt(A/π))
		vtkNew<vtkFloatArray> scale;
		scale->SetNumberOfComponents(1);
		scale->SetNumberOfTuples(N);
		scale->SetName("radius");
		for (int i = 0; i < N; ++i) {
			float r = std::sqrt(area_array[i] / pi);
			scale->SetValue(i, r);
		}

		// Set up source geometry: a unit disk
		vtkNew<vtkDiskSource> disk;
		disk->SetInnerRadius(0.0);     // solid disk
		disk->SetOuterRadius(1.0);     // will be scaled by Scale array
		disk->SetRadialResolution(1);  // low res for performance
		disk->SetCircumferentialResolution(20); // higher value for better visuals
		vtkNew<vtkTransform> transform;
		transform->RotateY(90);
		vtkNew<vtkTransformPolyDataFilter> transformFilter;
		transformFilter->SetInputConnection(disk->GetOutputPort());
		transformFilter->SetTransform(transform);
		transformFilter->Update();

		// Set up glyphs
		vtkNew<vtkPolyData> glyph_input;
		glyph_input->SetPoints(points);
		// Set vectors explicitly
		glyph_input->GetPointData()->AddArray(normals);
		glyph_input->GetPointData()->SetVectors(normals);  // <- Add this!
		glyph_input->GetPointData()->SetNormals(normals);
		glyph_input->GetPointData()->AddArray(scale);

		//Set height
		vtkNew<vtkFloatArray> height;
		height->SetNumberOfComponents(1);
		height->SetNumberOfTuples(N);
		height->SetName("height");
		for (int i = 0; i < N; ++i) {
			height->SetValue(i, height_array[i]);
		}
		glyph_input->GetPointData()->AddArray(height);

		vtkNew<vtkGlyph3D> glyphs;
		//glyphs->SetSourceConnection(disk->GetOutputPort());
		glyphs->SetSourceConnection(transformFilter->GetOutputPort());
		glyphs->SetInputData(glyph_input);
		glyphs->SetVectorModeToUseVector();       // <- Use "Normals" array
		glyphs->SetScaleModeToScaleByScalar();
		glyphs->SetColorModeToColorByScalar();
		glyphs->SetScaleFactor(1.0);
		glyphs->OrientOn();
		glyphs->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "radius");
		glyphs->Update();

		// Write to .vtu
		vtkNew<vtkXMLPolyDataWriter> writer;
		writer->SetFileName(path.string().c_str());
		writer->SetInputData(glyphs->GetOutput());
		writer->SetDataModeToBinary();
		writer->Write();
	}

    void OutputPointCloudDataAsVTP(const Array<Vector3>& positions, const std::vector<std::pair<Array<real>, std::string>> scalar_fields, const std::vector<std::pair<Array<Vector3>, std::string>> vector_fields, const fs::path& path) {
        Info("Outputting point cloud to VTP file: {}", path.string());

        vtkNew<vtkPolyData> poly_data;
        vtkNew<vtkPoints> vtk_points;
        vtk_points->SetNumberOfPoints(positions.size());

        for (vtkIdType i = 0; i < positions.size(); ++i) {
            const auto& p = positions[i];
            vtk_points->SetPoint(i, p[0], p[1], p[2]);
        }
        poly_data->SetPoints(vtk_points);

        // Add scalar fields
        for (const auto& [data, name] : scalar_fields) {
            if (data.size() != positions.size()) {
                Error("Scalar field '{}' size mismatch with positions\n", name);
                continue;
            }
            vtkNew<vtkFloatArray> arr;
            arr->SetName(name.c_str());
            arr->SetNumberOfComponents(1);
            arr->SetNumberOfTuples(positions.size());
            for (vtkIdType i = 0; i < positions.size(); ++i) {
                arr->SetValue(i, data[i]);
            }
            poly_data->GetPointData()->AddArray(arr);
        }

        // Add vector fields
        for (const auto& [data, name] : vector_fields) {
            if (data.size() != positions.size()) {
                Error("Vector field '{}' size mismatch with positions\n", name);
                continue;
            }
            vtkNew<vtkFloatArray> arr;
            arr->SetName(name.c_str());
            arr->SetNumberOfComponents(3);
            arr->SetNumberOfTuples(positions.size());
            for (vtkIdType i = 0; i < positions.size(); ++i) {
                const auto& v = data[i];
                arr->SetTuple3(i, v[0], v[1], v[2]);
            }
            poly_data->GetPointData()->AddArray(arr);
        }

        // Write to file
        vtkNew<vtkXMLPolyDataWriter> writer;
        writer->SetFileName(path.string().c_str());
        writer->SetInputData(poly_data);
        writer->SetDataModeToBinary();  // smaller file size
        writer->Write();
    }
}