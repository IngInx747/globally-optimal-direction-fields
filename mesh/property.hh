#ifndef OPENMESH_PROPERTY_EXTENSION_HH
#define OPENMESH_PROPERTY_EXTENSION_HH

#include <OpenMesh/Core/Utils/PropertyManager.hh>

namespace OpenMesh
{
/// 
/// @relates ConstPropertyViewer
/// @pre Property with the name \p propname of matching type exists.
/// @throws std::runtime_error if no property with the name \p propname of matching type exists.
/// @param mesh The mesh on which the property is created
/// @param propname The name of the created property
/// @tparam ElementT Element type of the created property, e.g. VertexHandle, HalfedgeHandle, etc.
/// @tparam T Value type of the created property, e.g., \p double, \p int, etc.
/// @returns A ConstPropertyViewer wrapping the property
/// 
template<typename ElementT, typename T>
inline ConstPropertyViewer<typename HandleToPropHandle<ElementT, T>::type>
getProperty(const PolyConnectivity &mesh, const char *propname)
{
    typename HandleToPropHandle<ElementT, T>::type prop;

    if (!mesh.get_property_handle(prop, propname))
    {
        std::ostringstream oss;
        oss << "Requested property handle \"" << propname << "\" does not exist.";
        throw std::runtime_error(oss.str());
    }

    return ConstPropertyViewer<decltype(prop)>(mesh, prop);
}

template<typename HandleT, typename ValueT>
inline const ValueT *getPropertyPtr(const PolyConnectivity &mesh, const char *propname)
{
    typename HandleToPropHandle<HandleT, ValueT>::type prop;
    if (!mesh.get_property_handle(prop, propname)) return nullptr;
    return mesh.property(prop).data();
}

template<typename HandleT, typename ValueT>
inline void removeProperty(PolyConnectivity &mesh, const char *propname)
{
    typename HandleToPropHandle<HandleT, ValueT>::type prop;
    if (mesh.get_property_handle(prop, propname)) mesh.remove_property(prop);
}

} // namespace OpenMesh

#endif