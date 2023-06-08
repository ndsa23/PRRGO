from django.urls import path
from . import views
from django.views.generic.base import TemplateView
urlpatterns = [
    path('go-viz/', TemplateView.as_view(template_name='go_viz.html')),
    path('submit/', views.cyto_json),
    path('download/', views.download)
]
