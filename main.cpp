#include "distance.h"
#include "timer.h"
#include "mst.cpp"
#include<sstream>
#include<unistd.h>
#include <sys/resource.h>
#include <sys/times.h>
#include <iostream>
#include <sstream>
#include<climits>
#include<iomanip>
#include "geodesic_algorithm_subdivision.h"
#define qtimes 100
//-------------------------------------
// MAIN
//-------------------------------------
int query_type=0;
int clpsize=6;
int algo_type=0;
int k=3;
int i;
int x[upper_poi];
FILE *fp;
int num_poi;
std::vector<int> poilist;
char prefix[255];
std::ofstream pathresult("path.txt", std::ios::out);
#ifndef WIN32
    double Time_preprocess=0;
    double  Time_dquery=0, Time_knnquery, Time_clpquery;
    double  Space_preprocess=0;
    double  Space_query=0;
    double errorbound_dis, errorbound_knn=0;
    struct rusage myTime_program_start, myTime_preprocess_end, myTime_query_begin, myTime_query_end;
#endif


double shortestpath_GB(int x, int y, geodesic::GeodesicAlgorithmExact& shortestpath){
    geodesic::SurfacePoint source(&mesh.vertices()[x]);
    geodesic::SurfacePoint dest(&mesh.vertices()[y]);
    std::vector<geodesic::SurfacePoint> sources;
    std::vector<geodesic::SurfacePoint> dests;
    sources.clear();
    dests.clear();
    sources.push_back(source);
    dests.push_back(dest);
    shortestpath.propagate_GB(sources, &dests);
    double dist;
    shortestpath.best_source(dest, dist);
    std::vector<geodesic::SurfacePoint> path;
    shortestpath.trace_back(dest, path);
    print_info_about_path(path);
    pathresult << std::endl << "Path from Vertex " << x << " to Vertex " << y << ":" << std::endl;
    for(int i=0;i<path.size();i++){
	    pathresult << path[i].getx() << " " << path[i].gety() << " " << path[i].getz() << std::endl;
    }
    return dist;
   }
void algo(geodesic::GeodesicAlgorithmExact &algorithm){
        std::vector<GeoNode*> AllPOI;
        AllPOI.clear();
        std::vector<std::pair<int, GeoNode*> > POIs;
        POIs.clear(); 
 //       fp=fopen("POINT.C","r");
        for(i=0;i<num_poi;i++)
        {
     //       if(fscanf(fp,"%d",&x[i]) != EOF)
            {
                GeoNode *n=new GeoNode(poilist[i],0);
                AllPOI.push_back(n);
                std::pair<int , GeoNode*> m(poilist[i], n);
                POIs.push_back(m);
 //               printf("x[%d]: %d\n",i,x[i]);
            }
        }
//       fclose(fp);
       printf("file read finished\n");
#ifndef WIN32
       getrusage(RUSAGE_SELF,&myTime_program_start);
#endif
        double radius=0;
        double distance;
        stx::btree<int, GeoNode*> bplusgeotree(POIs.begin(), POIs.end());
	    geodesic::SurfacePoint source(&mesh.vertices()[poilist[0]]);
        std::vector<geodesic::SurfacePoint> all_sources(1,source);
	    double const distance_limit = 0;
        geodesic::SurfacePoint target(&mesh.vertices()[poilist[1]]);
	    std::vector<geodesic::SurfacePoint> stop_points(1, target);
        for(i=2;i<num_poi;i++){
            stop_points.push_back(geodesic::SurfacePoint(&mesh.vertices()[poilist[i]]));
        }
	    algorithm.propagate(all_sources, distance_limit, &stop_points);
        for(i=1;i<num_poi;i++){
            geodesic::SurfacePoint p(&mesh.vertices()[poilist[i]]);
            algorithm.best_source(p,distance);
            radius= std::max(distance, radius);
        }
        std::cout<<"radius stopdis: "<< radius <<" "<< algorithm.distance_stopped()<<std::endl;
        GeoNode rootGeo(poilist[0],radius);
        stop_points.clear();
        BuildGeoTree(rootGeo,poi,bplusgeotree,algorithm);
        //std::cout<<"Tree Completed."<<std::endl;
	
     //   calculate_size(rootGeo);
      //  PrintGeoTree(rootGeo);
        
        generate_pair_geo(rootGeo, rootGeo, algorithm);
        std::cout<<"Preprocessing Finished."<<std::endl;
       /* std::ostringstream cmd;
        cmd.clear();
        cmd.str("");
        cmd << "/proc/" << getpid() << "/statm";
        fp=fopen(cmd.str().c_str(),"r");
        fscanf(fp,"%f", Space_preprocess);
        fclose(fp);*/
       // Space_query=2*geopairs.size()*(sizeof(GeoNode));
        Space_query=geopairs.size()*(sizeof(float)+sizeof(int)+sizeof(bool)+2*sizeof(char))/(1024.0*1024.0);
    //    Space_preprocess=2*geopairs.size()*sizeof(GeoNode)+sizeof(mesh)+sizeof(algorithm);
#ifndef WIN32
       getrusage(RUSAGE_SELF,&myTime_preprocess_end);
       Time_preprocess=myTime_preprocess_end.ru_utime.tv_sec-myTime_program_start.ru_utime.tv_sec;
#endif
     //  geopairs.sort(compair_geopair);
     //  PrintGeoPair();
       std::ofstream spanner("spanner.txt",std::ios::out);
       spanner << mesh.vertices().size() << " " << geopairs.size() << std::endl;
       for(std::list<GeoPair *>::iterator iten=geopairs.begin();iten!=geopairs.end();iten++){
            geopairsvector.push_back(*iten);
	    spanner << (*iten)->node1->index << " " << (*iten)->node2->index << " " << int((*iten)->distance) << std::endl;
       } 
/*       for(int m=0;m<geopairsvector.size();m++){
           for(std::list<int>::iterator bite=geopairsvector[m]->node1->mcode.begin();bite!=geopairsvector[m]->node1->mcode.end();bite++)
               std::cout<<*bite;
           std::cout<<", ";
           for(std::list<int>::iterator bite=geopairsvector[m]->node2->mcode.begin();bite!=geopairsvector[m]->node2->mcode.end();bite++)
               std::cout<<*bite;
           std::cout<<std::endl;
       }*/

/*
 * KNN Query
 * */
/*
#ifndef WIN32
       getrusage(RUSAGE_SELF,&myTime_query_begin);
#endif
       distance_return=knn(*geonodevector[qindex1],k);
       std::cout<<"distance_return:"<<distance_return<<std::endl;
#ifndef WIN32
       getrusage(RUSAGE_SELF,&myTime_query_end);
       Time_knnquery+=myTime_query_end.ru_utime.tv_usec-myTime_query_begin.ru_utime.tv_usec;
#endif
       realdistance=surface_knn(geonodevector[qindex1]->index,k);
       std::cout<<"realdistance:"<<realdistance<<std::endl;
       if(realdistance!=0&&realdistance<distance_return)errorbound_knn+=std::abs((double)((distance_return/realdistance)-(double)1.0));
       std::cout<<errorbound_knn<<std::endl;
*/
       /*
        * B-closest pair
        * */
/*
#ifndef WIN32
       getrusage(RUSAGE_SELF,&myTime_query_begin);
#endif
       std::vector<GeoNode*> v1;
       v1.clear();
       for(i=0;i<clpsize/2;i++){
           v1.push_back(geonodevector[qindex1]);
           qindex1=(qindex1+1)%geonodevector.size();
       }
       std::vector<GeoNode*> v2;
       v2.clear();
       for(i=0;i<clpsize/2;i++){
           v2.push_back(geonodevector[qindex2]);
           qindex2=(qindex2+1)%geonodevector.size();
       }*/
   //    k_closest_pairs(v1,v2,k);
//#ifndef WIN32
//       getrusage(RUSAGE_SELF,&myTime_query_end);
//       Time_clpquery+=myTime_query_end.ru_utime.tv_usec-myTime_query_begin.ru_utime.tv_usec;
//#endif
//       errorbound_knn/=qtimes;

        DeleteGeoTree(rootGeo);
}
/*-------------------------------------------------------------------------
 * ------------------------MAIN--------------------------------------------
 *  -------------------------------------------------------------------------
 *
 * */


int main(int argc, char **argv) 
{
	if(argc < 2)
	{
		std::cout << "usage: mesh_file_name " << std::endl; //try: "hedgehog_mesh.txt 3 14" or "flat_triangular_mesh.txt 1"
		return 0;
	}

    s = atof(argv[2]);
	bool success = geodesic::read_mesh_from_file(argv[1],points,faces);
	if(!success)
	{
		std::cout << "something is wrong with the input file" << std::endl;
		return 0;
	}

    strcpy(prefix, argv[1]);
	mesh.initialize_mesh_data(points, faces);		//create internal mesh data structure including edges

    geodesic::GeodesicAlgorithmExact algorithm(&mesh);	//create exact algorithm for the mesh
    std::cout << s << std::endl;

	//geodesic::GeodesicAlgorithmSubdivision subdivision_algorithm(&mesh,2);	//with subdivision_level=0 this algorithm becomes Dijkstra, with subdivision_level->infinity it becomes exact
    // WRITE THE RANDOM POINT FILE 
/*    fp = fopen("POINT.C","w");
    if ( fp == NULL )
    {
        puts ( "Cannot open file" );
        exit(1);
    }*/
    //std::ofstream out("poilist.txt", std::ios::out);
    std::ifstream input("poilist.txt", std::ios::in);
    input >> num_poi;
    poilist.resize(num_poi);
    for(int i=0;i<num_poi;i++) input >> poilist[i];
//    x[0]=randn(mesh.vertices().size());
   // out << poi << std::endl;
  //  out << x[0] << " ";
/*    for(i=0;i<poi;i++)
    {
        bool key=true;
        while(key){
            key=false;
            x[i]=randn(mesh.vertices().size());
            out << x[i] << " ";
            for(int j=0;j<i;j++){
                if(x[i]==x[j]){
                    key=true;
                    break;
                }
            }
        }
 //       fprintf(fp,"%d\n",x[i]);
    }
    out.close();*/
 //   fclose(fp);
// READ THE RANDOM POINT FILE AND ASSIGN TO ROOT Node
//    fp=fopen("output.txt","w");

//    fclose(fp);
    std::cout << "Finished reading POI list." << std::endl;
    pairs=0;
    pairvector.clear();
    nodevector.clear();
    geopairs.clear();
    geopairsvector.clear();
    geonodevector.clear();
    graphnodevector.clear();

    algo(algorithm);
    mst_preorder(poilist[0], poilist[num_poi-1]);

    ifstream orderfile("order.txt");
    int src, des;
    orderfile >> src;
    for(int i = 0;i<num_poi-1;i++){
	orderfile >> des;
	std::cout << src << " " << des << " ";
        shortestpath_GB(src, des, algorithm);
	src = des;
    }
}
