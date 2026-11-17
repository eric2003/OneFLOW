#include <vtkCubeSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>

int main(int argc, char *argv[])
{
	vtkNew<vtkCubeSource> cube;

	vtkNew<vtkPolyDataMapper> cubeMapper;
	cubeMapper->SetInputConnection( cube->GetOutputPort() );

	vtkNew<vtkActor> cubeActor;
	cubeActor->SetMapper( cubeMapper );

	cubeActor->GetProperty()->SetColor(1.0, 0.0, 0.0 );

	vtkNew<vtkRenderer> renderer;
	renderer->SetBackground( 0.0, 0.0, 1.0 );
	renderer->AddActor(cubeActor);

	vtkNew<vtkRenderWindow> renderWindow;
	renderWindow->AddRenderer(renderer);
	renderWindow->SetSize(400, 400);

	vtkNew<vtkRenderWindowInteractor> interactor;
	interactor->SetRenderWindow( renderWindow );
	interactor->Initialize();

	renderWindow->SetWindowName("OneFLOW CFD VTK Cube Test");
	renderWindow->Render();
	interactor->Start();

	return 0;
}
