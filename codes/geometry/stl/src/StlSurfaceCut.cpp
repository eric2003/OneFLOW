/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "StlSurfaceCut.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <iterator>
#include <stdexcept>
#include <utility>
#include <vector>

namespace ONEFLOW {
namespace {

using Point = std::array<double, 3>;
using Point2 = std::array<double, 2>;
using Edge = std::pair<std::size_t, std::size_t>;

Point Add(const Point& a, const Point& b) { return {a[0]+b[0],a[1]+b[1],a[2]+b[2]}; }
Point Sub(const Point& a, const Point& b) { return {a[0]-b[0],a[1]-b[1],a[2]-b[2]}; }
Point Mul(const Point& a, double s) { return {a[0]*s,a[1]*s,a[2]*s}; }
double Dot(const Point& a, const Point& b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
Point Cross(const Point& a, const Point& b)
{
    return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};
}
double Norm(const Point& a) { return std::sqrt(Dot(a,a)); }
Point Unit(const Point& a, const char* label)
{
    const double length = Norm(a);
    if (!std::isfinite(length) || length <= 1.0e-15)
        throw std::invalid_argument(std::string("Invalid ") + label + " vector");
    return Mul(a, 1.0/length);
}
double Distance2(const Point& a, const Point& b)
{
    const Point d = Sub(a,b);
    return Dot(d,d);
}
double SignedDistance(const Point& p, const StlCutPlane& plane, const Point& normal)
{
    return Dot(Sub(p,plane.point),normal);
}

struct WeldKey
{
    long long x,y,z;
    bool operator<(const WeldKey& other) const
    {
        return std::tie(x,y,z) < std::tie(other.x,other.y,other.z);
    }
};

WeldKey Key(const Point& p, double tolerance)
{
    return {
        static_cast<long long>(std::floor(p[0]/tolerance)),
        static_cast<long long>(std::floor(p[1]/tolerance)),
        static_cast<long long>(std::floor(p[2]/tolerance))
    };
}

struct WeldedPoints
{
    std::vector<Point> points;
    std::map<WeldKey,std::vector<std::size_t>> buckets;
    std::size_t Insert(const Point& point, double tolerance)
    {
        const WeldKey key = Key(point,tolerance);
        std::size_t best = points.size();
        for (long long x=-1;x<=1;++x)
        for (long long y=-1;y<=1;++y)
        for (long long z=-1;z<=1;++z)
        {
            const auto found=buckets.find({key.x+x,key.y+y,key.z+z});
            if(found==buckets.end()) continue;
            for(std::size_t id:found->second)
                if(Distance2(point,points[id])<tolerance*tolerance) best=std::min(best,id);
        }
        if(best==points.size())
        {
            best=points.size();
            points.push_back(point);
            buckets[key].push_back(best);
        }
        return best;
    }
};

std::vector<Point> ClipPolygon(
    const std::vector<Point>& polygon, const StlCutPlane& plane,
    const Point& normal, double tolerance)
{
    std::vector<Point> output;
    if(polygon.empty()) return output;
    auto distance=[&](const Point& p)
    {
        const double d=SignedDistance(p,plane,normal);
        return plane.keepPositive ? d : -d;
    };
    Point previous=polygon.back();
    double previousDistance=distance(previous);
    bool previousInside=previousDistance>=-tolerance;
    for(const Point& current:polygon)
    {
        const double currentDistance=distance(current);
        const bool currentInside=currentDistance>=-tolerance;
        if(currentInside!=previousInside)
        {
            const double denominator=previousDistance-currentDistance;
            if(std::abs(denominator)>std::numeric_limits<double>::epsilon())
            {
                const double t=std::clamp(previousDistance/denominator,0.0,1.0);
                output.push_back(Add(previous,Mul(Sub(current,previous),t)));
            }
        }
        if(currentInside) output.push_back(current);
        previous=current;
        previousDistance=currentDistance;
        previousInside=currentInside;
    }
    std::vector<Point> clean;
    for(const Point& p:output)
        if(clean.empty() || Distance2(clean.back(),p)>tolerance*tolerance*0.01)
            clean.push_back(p);
    if(clean.size()>1 && Distance2(clean.front(),clean.back())<=tolerance*tolerance*0.01)
        clean.pop_back();
    return clean;
}

void AddPolygon(StlMesh& mesh,const std::vector<Point>& polygon)
{
    if(polygon.size()<3) return;
    for(std::size_t i=1;i+1<polygon.size();++i)
    {
        StlTriangle tri;
        tri.vertices={polygon[0],polygon[i],polygon[i+1]};
        const Point n=Cross(Sub(tri.vertices[1],tri.vertices[0]),Sub(tri.vertices[2],tri.vertices[0]));
        const double length=Norm(n);
        if(length>0.0) tri.normal=Mul(n,1.0/length);
        mesh.triangles.push_back(tri);
    }
}

std::vector<Point> ExtractLoop(
    const std::vector<Edge>& edges,const std::vector<Point>& points,const char* label)
{
    if(edges.size()<3) throw std::runtime_error(std::string(label)+" cut did not produce a closed loop");
    std::map<std::size_t,std::set<std::size_t>> adjacency;
    for(const Edge& e:edges){adjacency[e.first].insert(e.second);adjacency[e.second].insert(e.first);}
    for(const auto& item:adjacency)
        if(item.second.size()!=2)
            throw std::runtime_error(std::string(label)+" cut is open, branched, or self-touching");
    const std::size_t start=adjacency.begin()->first;
    std::size_t previous=std::numeric_limits<std::size_t>::max();
    std::size_t current=start;
    std::vector<std::size_t> ids;
    do
    {
        ids.push_back(current);
        const auto& neighbors=adjacency.at(current);
        auto it=neighbors.begin();
        std::size_t next=*it;
        if(next==previous) next=*std::next(it);
        previous=current;
        current=next;
        if(ids.size()>adjacency.size())
            throw std::runtime_error(std::string(label)+" cut loop traversal did not close");
    }while(current!=start);
    if(ids.size()!=adjacency.size())
        throw std::runtime_error(std::string(label)+" cut produced multiple loops");
    std::vector<Point> loop;
    for(std::size_t id:ids) loop.push_back(points[id]);
    return loop;
}

double Orient2(const Point2& a,const Point2& b,const Point2& c)
{
    return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]);
}

std::vector<std::array<Point,3>> TriangulateCap(
    const std::vector<Point>& loop,const Point& planeNormal,
    const Point& outward,double tolerance)
{
    Point reference{1,0,0};
    if(std::abs(Dot(reference,planeNormal))>0.8) reference={0,1,0};
    const Point u=Unit(Cross(reference,planeNormal),"cap projection basis");
    const Point v=Cross(planeNormal,u);
    std::vector<Point2> projected;
    projected.reserve(loop.size());
    for(const Point& p:loop)
    {
        const Point q=Sub(p,loop.front());
        projected.push_back({Dot(q,u),Dot(q,v)});
    }
    double area=0.0;
    for(std::size_t i=0;i<projected.size();++i)
    {
        const Point2& a=projected[i];
        const Point2& b=projected[(i+1)%projected.size()];
        area+=a[0]*b[1]-b[0]*a[1];
    }
    if(std::abs(area)<=tolerance*tolerance)
        throw std::runtime_error("Cut boundary encloses zero area");
    const double winding=area>0.0?1.0:-1.0;
    std::vector<std::size_t> remaining(loop.size());
    for(std::size_t i=0;i<remaining.size();++i) remaining[i]=i;
    std::vector<std::array<Point,3>> triangles;
    auto inside=[&](const Point2& p,const Point2& a,const Point2& b,const Point2& c)
    {
        const double eps=tolerance*tolerance;
        return winding*Orient2(a,b,p)>=-eps &&
               winding*Orient2(b,c,p)>=-eps &&
               winding*Orient2(c,a,p)>=-eps;
    };
    while(remaining.size()>3)
    {
        bool found=false;
        for(std::size_t pos=0;pos<remaining.size();++pos)
        {
            const std::size_t a=remaining[(pos+remaining.size()-1)%remaining.size()];
            const std::size_t b=remaining[pos];
            const std::size_t c=remaining[(pos+1)%remaining.size()];
            if(winding*Orient2(projected[a],projected[b],projected[c])<=tolerance*tolerance)
                continue;
            bool contains=false;
            for(std::size_t id:remaining)
            {
                if(id==a||id==b||id==c) continue;
                if(inside(projected[id],projected[a],projected[b],projected[c]))
                {contains=true;break;}
            }
            if(contains) continue;
            triangles.push_back({loop[a],loop[b],loop[c]});
            remaining.erase(remaining.begin()+static_cast<std::ptrdiff_t>(pos));
            found=true;
            break;
        }
        if(!found) throw std::runtime_error("Cut loop cannot be triangulated as a simple polygon");
    }
    triangles.push_back({loop[remaining[0]],loop[remaining[1]],loop[remaining[2]]});
    for(auto& tri:triangles)
    {
        const Point normal=Cross(Sub(tri[1],tri[0]),Sub(tri[2],tri[0]));
        if(Dot(normal,outward)<0.0) std::swap(tri[1],tri[2]);
    }
    return triangles;
}

void AppendTriangles(StlMesh& target,const StlMesh& source)
{
    target.triangles.insert(target.triangles.end(),source.triangles.begin(),source.triangles.end());
}

} // namespace

StlSurfaceCutResult CutStlSurfaceBetweenPlanes(
    const StlMesh& source,const StlSurfaceCutOptions& options)
{
    if(!std::isfinite(options.tolerance)||options.tolerance<=0.0)
        throw std::invalid_argument("STL cutting tolerance must be positive and finite");
    const Point inletNormal=Unit(options.inlet.normal,"inlet plane normal");
    const Point outletNormal=Unit(options.outlet.normal,"outlet plane normal");
    const auto inputTopology=AnalyzeStlSurface(source,options.tolerance);
    if(!inputTopology.IsClosedTwoManifold()||inputTopology.connectedComponentCount!=1||
       inputTopology.signedVolume<=options.tolerance*options.tolerance*options.tolerance)
        throw std::invalid_argument("STL cutting requires one closed, connected, outward-oriented manifold shell");

    StlSurfaceCutResult result;
    result.wall.header="oneflow_stl_cut_wall";
    result.inletCap.header="oneflow_stl_cut_inlet";
    result.outletCap.header="oneflow_stl_cut_outlet";
    result.closedSurface.header="oneflow_stl_cut_closed_surface";
    for(const StlTriangle& triangle:source.triangles)
    {
        std::vector<Point> polygon(triangle.vertices.begin(),triangle.vertices.end());
        polygon=ClipPolygon(polygon,options.inlet,inletNormal,options.tolerance);
        polygon=ClipPolygon(polygon,options.outlet,outletNormal,options.tolerance);
        AddPolygon(result.wall,polygon);
    }
    if(result.wall.triangles.empty())
        throw std::runtime_error("Cutting planes removed the complete STL surface");

    WeldedPoints welded;
    std::map<Edge,std::size_t> edgeUses;
    for(const StlTriangle& tri:result.wall.triangles)
    {
        std::array<std::size_t,3> ids{};
        for(std::size_t i=0;i<3;++i) ids[i]=welded.Insert(tri.vertices[i],options.tolerance);
        for(std::size_t i=0;i<3;++i)
        {
            const std::size_t a=ids[i],b=ids[(i+1)%3];
            ++edgeUses[std::minmax(a,b)];
        }
    }
    std::vector<Edge> inletEdges,outletEdges;
    for(const auto& [edge,count]:edgeUses)
    {
        if(count!=1) continue;
        const Point midpoint=Mul(Add(welded.points[edge.first],welded.points[edge.second]),0.5);
        if(std::abs(SignedDistance(midpoint,options.inlet,inletNormal))<=options.tolerance*4.0)
            inletEdges.push_back(edge);
        else if(std::abs(SignedDistance(midpoint,options.outlet,outletNormal))<=options.tolerance*4.0)
            outletEdges.push_back(edge);
        else
            throw std::runtime_error("Clipped STL contains an unclassified boundary edge");
    }
    result.inletLoop=ExtractLoop(inletEdges,welded.points,"Inlet");
    result.outletLoop=ExtractLoop(outletEdges,welded.points,"Outlet");

    if(options.capBoundaries)
    {
        const Point inletOutward=Mul(inletNormal,options.inlet.keepPositive? -1.0:1.0);
        const Point outletOutward=Mul(outletNormal,options.outlet.keepPositive? -1.0:1.0);
        for(const auto& tri:TriangulateCap(result.inletLoop,inletNormal,inletOutward,options.tolerance))
        {
            StlTriangle item; item.vertices=tri;
            const Point n=Unit(Cross(Sub(tri[1],tri[0]),Sub(tri[2],tri[0])),"inlet cap normal");
            item.normal=n; result.inletCap.triangles.push_back(item);
        }
        for(const auto& tri:TriangulateCap(result.outletLoop,outletNormal,outletOutward,options.tolerance))
        {
            StlTriangle item; item.vertices=tri;
            const Point n=Unit(Cross(Sub(tri[1],tri[0]),Sub(tri[2],tri[0])),"outlet cap normal");
            item.normal=n; result.outletCap.triangles.push_back(item);
        }
    }
    result.closedSurface=result.wall;
    AppendTriangles(result.closedSurface,result.inletCap);
    AppendTriangles(result.closedSurface,result.outletCap);
    if(options.capBoundaries)
    {
        result.topology=AnalyzeStlSurface(result.closedSurface,options.tolerance);
        if(!result.topology.IsClosedTwoManifold()||result.topology.connectedComponentCount!=1||
           result.topology.signedVolume<=options.tolerance*options.tolerance*options.tolerance)
            throw std::runtime_error("Capped STL failed closed-manifold or positive-volume validation");
    }
    return result;
}

} // namespace ONEFLOW
