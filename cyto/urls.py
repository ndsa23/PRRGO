from django.urls import path
from . import views
from django.views.generic.base import TemplateView
urlpatterns = [
    path('go-viz/', TemplateView.as_view(template_name='go_viz.html')),
    path('go-viz/<str:keyword>.json/<int:filter_n>/', views.cyto_json, name='cyto_json'),
]
