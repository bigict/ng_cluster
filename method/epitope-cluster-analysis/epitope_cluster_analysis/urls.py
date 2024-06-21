from cluster2.views import CSVResultView
from django.conf.urls import url
from django.views.generic import TemplateView
from cluster2 import views
urlpatterns =  [
    #url(r'^$', view='main_view', name='cluster2-main-view'),
    url(r'^$', views.get_name, name='cluster2'),
#    url(r'^result/$', view='result_view', name='result2-view'),
    url(r'^results', views.process, name='results-view'),
    url(r'^step_2$', views.add_session2, name='add_session2'),
    url(r'^step_3$', views.add_session3, name='add_session3'),
    url(r'^result$', views.result, name='result'),
    url(r'^step3', views.step3, name='step3-view'),
#    url(r'^results/$', TemplateView.as_view(template_name='cluster2/results.html'), name='cluster2-result'),
    url(r'^result_in_text/$', views.text_result_view, name='text_result-view'),
    url(r'^result_in_csv/$', CSVResultView.as_csv, {}, name='csv_result-view'),
    url(r'^visualization/$', TemplateView.as_view(template_name='cluster2/visualization.html'), name='cluster2-visualization'),
    url(r'^help/$', TemplateView.as_view(template_name='cluster2/help.html'), name='cluster2-help'),
    url(r'^example/$', TemplateView.as_view(template_name='cluster2/example.html'), name='cluster2-example'),
    url(r'^reference/$', TemplateView.as_view(template_name='cluster2/reference.html'), name='cluster2-reference'),
    url(r'^download/$', TemplateView.as_view(template_name='cluster2/download.html'), name='cluster2-download'),
    url(r'^contact/$', TemplateView.as_view(template_name='cluster2/contact.html'), name='cluster2-contact'),
]
