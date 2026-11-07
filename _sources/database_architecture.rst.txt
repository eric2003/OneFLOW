Database Architecture and Data Access Contracts
=================================================

Purpose
-------

The database subsystem provides two related but different runtime stores:

* **Named parameters** hold configuration and other scalar or array values under
  string keys. They are represented by :cpp:class:`DataBase`,
  :cpp:class:`DataPara`, :cpp:class:`DataEntry`, and
  :cpp:class:`DataObject`.
* **Field registrations** associate names with pointer wrappers for solver or
  mesh field memory. The registry is owned by :cpp:class:`DataBase`, but the
  field arrays referenced by it are managed by their respective owners.

This distinction matters: owning a pointer wrapper is not the same as owning
the numerical field storage to which that wrapper refers.

Ownership model
---------------

The current ownership chain for named parameters is:

.. code-block:: text

   DataBase
     +-- unique_ptr<DataPara>
           +-- unordered_map<string, unique_ptr<DataEntry>>
                 +-- unique_ptr<DataObject>
                       +-- TDataObject<T>
                             +-- vector<T>

Each level has one clear owner:

* `DataBase` owns the parameter container and field registry.
* `DataPara` owns each named `DataEntry`.
* `DataEntry` owns its type-erased `DataObject`.
* `TDataObject<T>` owns the actual `std::vector<T>` values.

The `unique_ptr` chain makes destruction deterministic. Removing an entry or
clearing `DataPara` destroys the corresponding entry and its stored values.
Pointers and references returned by lookup functions are borrowed; they must
not outlive the entry and must not be retained across removal or clearing.

The global database is itself held by a `std::unique_ptr<DataBase>`. The
global accessor pair has intentionally different contracts:

* `GetGlobalDataBase()` is nullable and is retained for compatibility with
  APIs that need to check initialization.
* `RequireGlobalDataBase()` expresses a required dependency and throws
  `std::logic_error` when the global instance is unavailable.

Parameter insertion and replacement
-----------------------------------

`ProcessData` converts parsed text into a concrete `TDataObject<T>`, wraps it
in a `DataEntry`, and transfers ownership to `DataPara::SetDataEntry`.

For a new key, `SetDataEntry` stores the incoming entry. For an existing key,
it calls `DataEntry::Copy`:

* The existing entry's name, type, and size remain unchanged.
* The incoming value is copied into the existing payload.
* A type or size mismatch throws `std::runtime_error`.
* The temporary incoming entry is destroyed automatically after the call.

This is an update-in-place value contract, not replacement of the entry object.
Code holding a borrowed pointer to an existing entry therefore does not see
that pointer invalidated by a successful same-key value update. Removal and
`Clear()`, however, do invalidate borrowed entry pointers.

The current text ingestion path in `ProcessData` explicitly constructs
`TDataObject<std::string>`, `TDataObject<int>`, and `TDataObject<Real>`.
The type constants declared elsewhere should not be interpreted as proof that
every possible type is supported by this ingestion path.

Type erasure and typed reads
----------------------------

`DataObject` is the polymorphic boundary between ownership/serialization code
and concrete value types. `TDataObject<T>` owns the vector and implements
copying, serialization, deserialization, and diagnostic dumping.

The low-level accessors are:

.. code-block:: cpp

   template < typename T >
   T GetDataValue( DataObject * dataObject, int iElement = 0 );

   template < typename T >
   T GetDataValue( const DataObject * dataObject, int iElement = 0 );

Both return a **value copy**. The const overload preserves read-only access
through the type-erased interface; it is not a cast-away-const convenience.
The default index selects the first element. The two overloads deliberately
have matching default arguments so mutable and const access have the same
call shape.

These are low-level accessors, not checked dynamic type queries. Their caller
must ensure all of the following:

* `dataObject` is non-null;
* `T` exactly matches the concrete `TDataObject<T>` used for storage;
* `iElement` is a valid zero-based index;
* the stored object contains at least one element when the default index is
  used.

The accessors use the type-erased `GetVoidPointer()` and a static cast. They
do not carry enough runtime type or length information to validate a caller's
template argument or index. Violating these preconditions can cause undefined
behavior. Prefer the higher-level `GetDataValue<T>(name, database)` overloads
when a named lookup is needed; those additionally check a null database and a
missing key, but still rely on the caller to request the correct stored type.

Const access contract
---------------------

The API provides const and mutable lookup at both levels:

* `DataBase*` lookup obtains a mutable `DataEntry` and `DataObject`.
* `const DataBase*` lookup obtains const entry/object references and calls
  the const `DataObject*` overload.
* `DataPara` provides matching const and mutable `FindDataEntry` overloads.

Constness prevents mutation through that access path. It does not make the
underlying global database immutable if other code retains a mutable handle.
The database is not documented here as thread-safe; concurrent mutation and
lookup require external synchronization unless a stronger contract is
established by the relevant runtime.

Field pointer registry
----------------------

`DataField` stores named `PointerWrap` objects. A wrapper provides access to
field memory managed elsewhere, such as a solver-owned vector or mesh array.
The database owns the registry entries and wrappers, but callers must preserve
the lifetime and address stability of the referenced field memory for as long
as the wrapper is used.

Use the typed field pointer/reference helpers instead of repeating raw
`PointerWrap` lookup logic. Mutable access should be used only where mutation
is intended; const storage should resolve to const pointers/references.

Serialization boundary
----------------------

Named entries serialize their name, type identifier, element count, and
payload through `DataBook`. The concrete `TDataObject<T>` handles element
serialization using the matching `HXWrite`/`HXRead` operations. This keeps
the database container independent of the concrete serialization of each
stored type.

Serialization support and text-ingestion support are separate concerns:
adding a type identifier alone does not establish that parsing, creation,
copying, and serialization are all implemented for that type.

Relationship to configuration ingestion
----------------------------------------

The configuration ingestion path composes the database as follows:

.. code-block:: text

   control.txt manifest
          |
          v
   ordered script paths
          |
          v
   production script parser
          |
          v
   ProcessData(name, text values, type, size)
          |
          v
   typed TDataObject<T> + DataEntry
          |
          v
   DataPara::SetDataEntry
          |
          v
   named runtime parameter lookup

The ingestion integration test verifies a real two-file manifest and parser
path, typed values, and later-file update of a repeated key. It exercises both
mutable and const scalar lookup. It does not prove the broader
`ReadControlInfo(caseDir)` orchestration, parallel broadcast, diagnostic
database dump, case setup/teardown, or sequential-case isolation. Those remain
separate runtime integration obligations tracked in
`doc/plans/setup-case-environment-integration-test.md`.

Maintenance rules
-----------------

When changing this subsystem:

#. Preserve the single-owner `unique_ptr` chain unless a measured use case
   requires shared ownership.
#. Document whether returned pointers/references are owning or borrowed and
   which operations invalidate them.
#. Keep mutable and const overloads behaviorally symmetric.
#. Treat type and element-count validation as explicit API responsibilities;
   do not assume a C++ cast validates the stored runtime type.
#. Keep parameter values separate from externally owned solver field arrays.
#. Add tests for observable contracts: missing keys, type/size mismatch,
   same-key updates, clear/removal invalidation boundaries, const reads, and
   serialization round trips.
#. Keep the configuration-ingestion tests distinct from full runtime lifecycle
   tests so failures identify the responsible boundary.

Source entry points
-------------------

* `codes/database/include/DataBase.h`
* `codes/database/include/DataPara.h`
* `codes/database/include/DataObject.h`
* `codes/database/src/DataBase.cpp`
* `codes/database/src/DataPara.cpp`
