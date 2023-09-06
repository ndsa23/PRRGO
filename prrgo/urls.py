from django.contrib import admin
from django.urls import include, path

urlpatterns = [
    path('', include(('cyto.urls', 'cyto'), namespace='cyto')),
    path('admin/', admin.site.urls),
]
