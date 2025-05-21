import React, { useState, useId, useDeferredValue, useTransition, memo } from "react";
import { cn } from "../../lib/utils";
import { ChevronDown, Check } from "lucide-react";

// Memo-ize the Select component to prevent unnecessary re-renders
const Select = memo(function Select({ 
  children, 
  defaultValue, 
  value, 
  onChange, 
  disabled, 
  placeholder, 
  ...props 
}) {
  const [isOpen, setIsOpen] = useState(false);
  const [selectedValue, setSelectedValue] = useState(value || defaultValue || "");
  const [selectedLabel, setSelectedLabel] = useState("");
  const [isPending, startTransition] = useTransition();
  const deferredSelectedValue = useDeferredValue(selectedValue);
  
  // Generate unique ID for this select component
  const id = useId();
  
  React.useEffect(() => {
    if (value !== undefined) {
      // Use transition for smoother UI updates when the value changes
      startTransition(() => {
        setSelectedValue(value);
      });
    }
  }, [value, startTransition]);
  
  // Extract children to find SelectItem components
  const items = React.Children.toArray(children)
    .filter(child => React.isValidElement(child) && child.type === SelectContent)
    .flatMap(content => 
      React.Children.toArray(content.props.children)
        .filter(child => React.isValidElement(child) && child.type === SelectItem)
    );
  
  React.useEffect(() => {
    // Find the label for the selected value
    const item = items.find(item => 
      React.isValidElement(item) && item.props.value === deferredSelectedValue
    );
    
    if (item && React.isValidElement(item)) {
      startTransition(() => {
        setSelectedLabel(item.props.children);
      });
    } else {
      startTransition(() => {
        setSelectedLabel("");
      });
    }
  }, [deferredSelectedValue, items, startTransition]);
  
  const handleValueChange = (newValue) => {
    if (onChange) {
      onChange(newValue);
    }
    
    startTransition(() => {
      setSelectedValue(newValue);
      setIsOpen(false);
    });
  };
  
  const handleOpen = () => {
    if (!disabled) {
      setIsOpen(prev => !prev);
    }
  };
  
  // Create context to share state with child components
  const selectContext = React.useMemo(() => ({
    open: isOpen,
    onOpenChange: setIsOpen,
    value: deferredSelectedValue,
    onValueChange: handleValueChange,
    disabled,
    id
  }), [isOpen, deferredSelectedValue, handleValueChange, disabled, id]);
  
  return (
    <SelectContext.Provider value={selectContext}>
      <div className={cn("relative", isPending && "cursor-wait opacity-80")} {...props}>
        {children}
      </div>
    </SelectContext.Provider>
  );
});

// Create context to share state between components
const SelectContext = React.createContext({});

const useSelectContext = () => React.useContext(SelectContext);

const SelectTrigger = React.forwardRef(({ className, children, ...props }, ref) => {
  const { open, onOpenChange, disabled, id } = useSelectContext();
  
  return (
    <button
      ref={ref}
      id={`${id}-trigger`}
      aria-haspopup="listbox"
      aria-expanded={open}
      aria-controls={`${id}-content`}
      disabled={disabled}
      onClick={() => onOpenChange(!open)}
      className={cn(
        "flex h-10 w-full items-center justify-between rounded-md border border-input bg-background px-3 py-2 text-sm ring-offset-background placeholder:text-muted-foreground focus:outline-none focus:ring-2 focus:ring-ring focus:ring-offset-2 disabled:cursor-not-allowed disabled:opacity-50",
        className
      )}
      {...props}
    >
      {children}
      <ChevronDown className={cn("h-4 w-4 opacity-50 transition-transform", open && "rotate-180")} />
    </button>
  );
});
SelectTrigger.displayName = "SelectTrigger";

const SelectValue = React.forwardRef(({ className, placeholder, ...props }, ref) => {
  const { value } = useSelectContext();
  
  return (
    <span 
      ref={ref} 
      className={cn("block truncate", className, !value && "text-muted-foreground")}
      {...props}
    >
      {props.children || placeholder}
    </span>
  );
});
SelectValue.displayName = "SelectValue";

const SelectContent = memo(function SelectContent({ className, children, ...props }) {
  const { open, id } = useSelectContext();
  
  if (!open) return null;
  
  return (
    <div
      id={`${id}-content`}
      role="listbox"
      className={cn(
        "absolute z-50 min-w-[8rem] w-full max-h-[--radix-select-content-available-height] mt-1 overflow-hidden rounded-md border bg-popover text-popover-foreground shadow-md animate-in fade-in-80",
        className
      )}
      {...props}
    >
      <div className="p-1 overflow-y-auto max-h-[300px]">
        {children}
      </div>
    </div>
  );
});

const SelectItem = React.forwardRef(({ className, children, value, ...props }, ref) => {
  const { value: selectedValue, onValueChange } = useSelectContext();
  const isSelected = selectedValue === value;
  
  return (
    <div
      ref={ref}
      role="option"
      aria-selected={isSelected}
      data-value={value}
      data-highlighted={isSelected ? true : undefined}
      className={cn(
        "relative flex w-full cursor-default select-none items-center rounded-sm py-1.5 pl-8 pr-2 text-sm outline-none focus:bg-accent focus:text-accent-foreground hover:bg-accent hover:text-accent-foreground data-[disabled]:pointer-events-none data-[disabled]:opacity-50",
        isSelected && "bg-accent text-accent-foreground",
        className
      )}
      onClick={() => onValueChange(value)}
      {...props}
    >
      <span className={cn(
        "absolute left-2 flex h-3.5 w-3.5 items-center justify-center", 
        !isSelected && "opacity-0"
      )}>
        <Check className="h-4 w-4" />
      </span>
      <span className="truncate">{children}</span>
    </div>
  );
});
SelectItem.displayName = "SelectItem";

export { Select, SelectTrigger, SelectValue, SelectContent, SelectItem };